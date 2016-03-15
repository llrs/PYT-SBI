#!/usr/bin/python3
# encoding: utf-8
"""This module provides a function to work with the online version of BLAST
provided by the NCBI.
@author: Leo, Lluís, Ferran"""

import logging

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq


def run_BLAST(query, db, output, size, filt=True):
    """Runs a psi-blast online.

    query is the file with the sequence
    db is the db to do the blast usually nr or pdb
    output is the file where the results will be stored
    size is the number of expected results.
    filt is the option to filter by genere,
    so only those sequence of different genere will be saved"""
    logging.info("Starting a blast from {} on {}.".format(query, db))
    record = SeqIO.read(query, format="fasta")
    result_handle = NCBIWWW.qblast('blastp', db, record.format("fasta"),
                                    service='psi',
                                    hitlist_size=size)

    blast_record = NCBIXML.read(result_handle)
    sq = len(record.seq)

    with open(output, 'w') as blast_record_out:
        SeqIO.write(record, blast_record_out, "fasta")
        genere_set = set()
#         for rounds in blast_record_out.rounds:
#             for alignment in rounds.alignments:
#                 print(alignment.title)
#                 for hsp in alignment.hsps:
#                     print(hsp.score)
        if filt:
            for alignment in blast_record.alignments:
                if "[" in alignment.title:
                    spiece = alignment.title.split('[')[-1].rstrip("]")
                    genere = spiece.split()[1]
                    for hsp in alignment.hsps:
                        percentage_identity = 100 * hsp.identities/sq
                        if percentage_identity > 30:
                            if genere not in genere_set:
                                blast_record_out.write(">{} {%.2f}\n{}\n".format(alignment.title, 
                                                                 percentage_identity, 
                                                                 hsp.sbjct))
                                genere_set.add(genere)
        else:
            for alignment in blast_record.alignments:
                print(alignment.title)
                for hsp in alignment.hsps:
                    SeqIO.write(Seq(hsp.sbjct), blast_record_out, format="fasta")


if __name__ == "__main__":
    run_BLAST('P01732.fasta', "pdb", "output_blast_2.out", 50, True)
