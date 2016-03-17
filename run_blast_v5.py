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

Entrez.email = "ferran.muinos@gmail.com"
Entrez.tool = "cozmic.py"

def run_BLAST(query, blast_type, db, size, filt=True):
    """Runs a blast online.

    query is the file with the sequence.
    db is the db to do the blast usually nr or pdb.
    output is the file where the results will be stored.
    size is the number of expected results.
    filt is the option to filter by genus,
    so only those sequence of different genus will be saved

    Returns a list of dictionaries with the id of the result of the blast"""
    logging.info("Starting a blast from {} on {} with {}.".format(
                                                        query, db, blast_type))
    if os.path.isfile(query):
        record = SeqIO.read(query, format="fasta")
        result_handle = NCBIWWW.qblast(blast_type, db, record.format("fasta"),
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
    """Generator downloading sequences from Entrez."""
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
