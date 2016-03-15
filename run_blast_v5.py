#!/usr/bin/python3
# encoding: utf-8
"""This module provides a function to work with the online version of BLAST
provided by the NCBI.
@author: Leo, Lluís, Ferran"""

import logging
import argparse

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
import os


def run_BLAST_original(query):
    blast_record_out = open('blast_record.fa','w')
    record = SeqIO.read(query, format="fasta")
    sz = len(record.seq)
    result_handle = NCBIWWW.qblast('blastp','nr', record.format("fasta"),
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
                percentage_identity = 100*hsp.identities/sz
                if percentage_identity > 30:
                    if genus not in genere_set:
#                         SeqIO.write(hsp.sbjct, blast_record_out, "fasta")
                        blast_record_out.write(">{}\n{}\n".format(alignment.title,
                                                                 hsp.sbjct))
                        genere_set.add(genus)


def run_BLAST(query, db, output, size, filt=True):
    """Runs a psi-blast online.

    query is the file with the sequence
    db is the db to do the blast usually nr or pdb
    output is the file where the results will be stored
    size is the number of expected results.
    filt is the option to filter by genere,
    so only those sequence of different genere will be saved"""
    logging.info("Starting a blast from {} on {}.".format(query, db))
    with open(output, 'w') as blast_record_out:
        if os.path.isfile(query):
            record = SeqIO.read(query, format="fasta")
            blast_record_out.write(">{}\n{}\n".format(record.id, record.seq))
#             SeqIO.write(record, blast_record_out, "fasta")
            result_handle = NCBIWWW.qblast('blastp', db, record.format("fasta"),
                                        service='psi',
                                        hitlist_size=size)
            sq = len(record.seq)
        else:
            result_handle = NCBIWWW.qblast('blastp', db, query,
                                        service='psi',
                                        hitlist_size=size)
        blast_record = NCBIXML.read(result_handle)

        genere_set = set()
        if filt:
            for alignment in blast_record.alignments:
                if "[" in alignment.title:
                    spiece = alignment.title.split('[')[-1].rstrip("]")
                    genere = spiece.split()[0]
                    for hsp in alignment.hsps:
                        percentage_identity = 100 * hsp.identities/sq
                        if percentage_identity > 30:
                            if genere not in genere_set:
#                                 print(genere)
                                blast_record_out.write(">{}\n{}\n".format(
                                alignment.title.split(">")[-1], hsp.sbjct))
#                                 SeqIO.write(hsp.sbjct, blast_record_out, "fasta")
                                genere_set.add(genere)
        else:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    blast_record_out.write(">{}\n{}\n".format(
                                alignment.title.split(">")[-1], hsp.sbjct))


if __name__ == "__main__":
    msg = 'Runs blast online.'
    argparser = argparse.ArgumentParser(description=msg,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    query = argparser.add_mutually_exclusive_group()
    query.add_argument("-i", help="Id of the sequence ")
    query.add_argument("-f", help="File with sequences")
    argparser.add_argument("-o", help="Output file")
    argparser.add_argument("-s", help="Set the number hits you want")
    argpaser.add_argument("-s", help="Filter by genus if possible",
                        action='store_true')
    run_BLAST('P01732', "pdb", "blast_record_5.out", 50, True)
#     run_BLAST_original('P01732.fasta')
#     out = run_online_blast("P01732", "pdb", 50)
#     handle_blast(out, "blast.out")
    
