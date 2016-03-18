#!/usr/bin/python3
# encoding: utf-8
"""This module provides a function to work with the online version of BLAST
provided by the NCBI.
@author: Leo, Lluís, Ferran"""

import logging
import argparse
import os
import urllib
import ftplib

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1
from Bio.Blast.Applications import NcbiblastpCommandline

import contact_map as cm
import plots

Entrez.email = "ferran.muinos@gmail.com"
Entrez.tool = "cozmic.py"


def run_BLAST(query, blast_type, db, size):
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
    return blast_record


def analyze_blast_result(blast_out, filt=True):
    """Classify the result of blast."""
    logging.info("Analysing the result of blast %s.", blast_out.query_id)
    sq = blast_out.query_length
    id_set = []
    genere_set = set()
    for alignment in blast_out.alignments:
        logging.debug("Analyzing alignment: %s", alignment.title)
        if filt:
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
                        msg_hsp = "Finishing high-scoring pair of an alignment"
                        logging.debug(msg_hsp)
        else:
            values = alignment.hit_id.split("|")
            id_set.append({values[0]: values[1], values[2]: values[3]})
    else:
        logging.debug("No more alignments left.")
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


def local_blast(query, blast_type, db, remote=True, **kwargs):
    """Function to run with the local blast program"""
    logging.info("Running blast locally with {} and {}".format(query, db))
    if remote:
        blast_cline = NcbiblastpCommandline(query=query, db=db,
                                            remote=True, out="blast.out",
                                            outfmt="5", evalue=0.001, **kwargs)
    else:
        blast_cline = NcbiblastpCommandline(query=query, db=db, outfmt="5",
                                            out="blast.out",
                                            evalue=0.001, **kwargs)
    print(blast_cline)
    stdout, stderr = blast_cline()
    logging.debug(stderr)
    blast_out = open("blast.out", "w")
    blast_record = NCBIXML.read(blast_out)
    return blast_record


if __name__ == "__main__":
    msg = 'Runs blast online.'
    args_helper = argparse.ArgumentDefaultsHelpFormatter
    argparser = argparse.ArgumentParser(description=msg,
                                        formatter_class=args_helper)

    argparser.add_argument("input",
                           help="Id of the sequence or file ")
    argparser.add_argument("output_file",
                           help="Output file")
    argparser.add_argument("type",
                           help="Type of blast to perform",
                           choices=["blastp", "blastn", "blast", "blastx",
                                    "tblastn", "tblastx"])
    choices_db = ["nr", "pdb", "swissprot", "refseq_protein", "pat",
                  "env_nr", "tsa_nr"]
    argparser.add_argument("db", 
                           help="Set the database to search on.",
                           choices=choices_db,
                           default="nr")
    argparser.add_argument("-s",
                           help="Set the number hits you want",
                           type=int,
                           default=200)
    argparser.add_argument("-f",
                           help="If present don't filter by genus",
                           action='store_false',
                           default=True)
    argparser.add_argument("-l",
                           help="""If present do a local blast on the db 
                           path""",
                           action="store_true", default=False)

    args = argparser.parse_args()

    blast_result = run_BLAST(args.input, args.type, args.db, args.s)
#     blast_result = local_blast(args.input, args.type, args.db)
    ides = analyze_blast_result(blast_result, args.f)
    ids = list(filter_ids(ides, "gi"))
    ides_pdb = list(filter_ids(ides, "pdb"))
    file_out = open(args.output_file, "w")

    for pdb in ides_pdb:
        try:
            pdb_file = plots.pdb_download(pdb, os.getcwd())
        except urllib.error.URLError:
            pass
        except ftplib.error_perm:
            pass
        else:
            parser = PDBParser(PERMISSIVE=1)
            structure = parser.get_structure(pdb, pdb_file)
            residues = cm.filter_residues(structure)
            residues_names = list(map(lambda x: x.get_id(), residues))
            seq = seq1(residues_names)
            seq.id = pdb
            SeqIO.write(seq, file_out, "fasta")
    else:
        if os.path.isfile(args.input):
            record = SeqIO.read(args.input, format="fasta")
            SeqIO.write(record, file_out, "fasta")
        else:
            ids.append(args.input)
        print(ids)
        SeqIO.write(retrive_sequence(ids), file_out, "fasta")
        file_out.close()
