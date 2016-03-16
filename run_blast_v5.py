#!/usr/bin/python3
# encoding: utf-8
"""This module provides a function to work with the online version of BLAST
provided by the NCBI.
@author: Leo, Lluís, Ferran"""

import logging
import argparse
import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "lluisrevillasa@gmail.com"
Entrez.tool = "correlatedmutations.py"


def run_BLAST_original(query):
    blast_record_out = open('blast_record.fa', 'w')
    record = SeqIO.read(query, format="fasta")
    sz = len(record.seq)
    result_handle = NCBIWWW.qblast('blastp', 'nr', record.format("fasta"),
                                service='psi',
                                hitlist_size=1000)
    blast_record = NCBIXML.read(result_handle)

#     blast_record_out.write(">{}\n{}\n".format(record.id, record.seq))
    SeqIO.write(record, blast_record_out, "fasta")
    genere_set = set()
    for alignment in blast_record.alignments:
        if '[' in alignment.title:
            spiece = alignment.title.split('[')[-1].rstrip("]")
            genus = spiece.split()[0]
            for hsp in alignment.hsps:
                percentage_identity = 100 * hsp.identities / sz
                if percentage_identity > 30:
                    if genus not in genere_set:
#                         SeqIO.write(hsp.sbjct, blast_record_out, "fasta")
                        blast_record_out.write(">{}\n{}\n".format(alignment.title,
                                                                 hsp.sbjct))
                        genere_set.add(genus)


def run_BLAST(query, blast_type, db, size, filt=True):
    """Runs a blast online.

    query is the file with the sequence.
    db is the db to do the blast usually nr or pdb.
    output is the file where the results will be stored.
    size is the number of expected results.
    filt is the option to filter by genus,
    so only those sequence of different genus will be saved"""
    logging.info("Starting a blast from {} on {}.".format(query, db))
    if os.path.isfile(query):
        record = SeqIO.read(query, format="fasta")
        result_handle = NCBIWWW.qblast(blast_type, db, record.seq,
                                    hitlist_size=size)
    else:
        result_handle = NCBIWWW.qblast(blast_type, db, query,
                                    hitlist_size=size)
    blast_record = NCBIXML.read(result_handle)

    sq = blast_record.query_length
    id_set = []
    genere_set = set()
    if filt:
        for alignment in blast_record.alignments:
            if "[" in alignment.title:
                spiece = alignment.title.split('[')[-1].rstrip("]")
                genere = spiece.split()[0]
                for hsp in alignment.hsps:
                    percentage_identity = 100 * hsp.identities / sq
                    if percentage_identity > 30:
                        if genere not in genere_set:
                            values = alignment.hit_id.split("|")
                            id_add = {values[0]: values[1],
                                      values[2]: values[3]}
                            id_set.append(id_add)
                            genere_set.add(genere)
    else:
        for alignment in blast_record.alignments:
            values = alignment.hit_id.split("|")
            id_set.append({values[0]: values[1], values[2]: values[3]})
    return id_set


def retrive_sequence(id_seqs):
    """Download sequence from Entrez"""
    logging.info("Downloading sequences from Entrez.")
    for id_seq in id_seqs:
        logging.debug("Downloading sequence {}.".format(id_seq))
        handle = Entrez.efetch(db="protein", id=id_seq, rettype="fasta",
                     retmode="text")
        yield SeqIO.read(handle, "fasta")


def filter_ids(ids, key):
    """Extract all the values of a shared key from a list."""
    logging.info("Extracting ids for %s.", key)
    return map(lambda x: x[key], ids)


if __name__ == "__main__":
    msg = 'Runs blast online.'
    argparser = argparse.ArgumentParser(description=msg,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    argparser.add_argument("i", help="Id of the sequence or file ")
    argparser.add_argument("output_file", help="Output file")
    argparser.add_argument("type", help="Type of blast to perform",
                           choices=["blastp", "blastn", "blast", "blastx",
                                    "tblastn", "tblastx"])
    choices_db = ["nr", "pdb", "swissprot", "refseq_protein", "pat",
                  "env_nr", "tsa_nr"]
    argparser.add_argument("db", help="Set the database to search on.",
                           choices=choices_db, default="nr")
    argparser.add_argument("-s", help="Set the number hits you want",
                        const=50, action='store_const')
    argparser.add_argument("-f", help="If present don't filter by genus",
                        action='store_false', default=True)

    args = argparser.parse_args()

    ides = run_BLAST(args.i, args.type, args.db, args.s, args.f)

    ids = list(filter_ids(ides, "gi"))
    file_out = open(args.output_file, "w")
    if os.path.isfile(args.i):
        record = SeqIO.read(args.i, format="fasta")
        SeqIO.write(record, file_out, "fasta")
    else:
        ids.append(args.i)
    SeqIO.write(retrive_sequence(ids), file_out, "fasta")
    file_out.close()
