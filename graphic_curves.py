#!/usr/bin/python3
# encoding: utf-8
'''
Provides the procedure to plot a graph with cutoffs against predicted
number of CM hits against precision to predict contact
Created on Mar 14, 2016
@author: Leo, Ferran, Lluís
'''
# Standard modules
import argparse
from distutils.spawn import find_executable
import matplotlib.pyplot as plt
import numpy as np

# Biopython modules
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq

# Non-standard modules
import contact_map as cm
import mutual_information as mut
import msa_caller as msa


def plot_twin_curves(cutoff_list, hit_list, precision_list):
    """Plots the cutoff_list x hit_list and cutoff_list x precision_list
    curves in the same plot with two different yaxis labels."""
    fig, ax1 = plt.subplots()
    ax1.plot(cutoff_list, hit_list, 'ob-')
    ax1.set_xlabel('Cutoff level L')
    # Make the y-axis label match the line color.
    ax1.set_ylabel('Predicted CM residue pairs', color='b')
    # Make the second plot with another axis
    ax2 = ax1.twinx()
    ax2.plot(cutoff_list, precision_list, 'or-')
    ax2.set_ylabel('Precision', color='r')
    plt.show()


def matrix_hits(binary_matrix):
    """Number of True cells above the diagonal in a binary matrix"""
    count = 0
    (n1, n2) = binary_matrix.shape
    for i in range(n1):
        for j in range(n2):
            if j > i:
                count += binary_matrix[i, j]
    return count


if __name__ == '__main__':
    # Argument parsing
    msg = 'Graphs with cutoff against CM hits and precision'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="protein structure file PDB")
    argparser.add_argument("-a",
                           help="""Atom to calculate distance with
                           CA: Carbon Alpha, CB: Carbon Beta. By default
                           is the minimum distance between residues""",
                           default="min", choices=["CA", "CB"])
    argparser.add_argument("-b",
                           help="Base of the logarithms.",
                           default=20)
    argparser.add_argument("-low",
                           help="Min entropy threshold.",
                           default=0.3)
    argparser.add_argument("-high",
                           help="Max entropy threshold.",
                           default=0.9)
    argparser.add_argument("-m",
                           help="Method for MSA.",
                           default="clustalw",
                           choices=["muscle", "t_coffee"])
    args = argparser.parse_args()

    # get PDB and distance matrix
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("test", args.file)
    residues = cm.filter_residues(structure)
    dist_matrix = cm.calc_dist_matrix(residues, args.a)
    cont_matrix = cm.contact_map(dist_matrix, args.a)
    # get the residue sequence from the PDB
    ppb = PPBuilder()
    seq = Seq('', generic_protein)
    for pp in ppb.build_peptides(structure):
        seq += pp.get_sequence()
    # keep sequence in a FASTA to call BLAST afterwards
    fd = open('temp.fa', 'w')
    fd.write(">target\n%s" % seq)
    fd.close()
    # call BLAST search: as of today the output is named "blast_record.fa"
    blst.run_BLAST('temp.fa')
    # get the MSA from the blast_record
    msa.call_msa_method(args.m, "blast_record.fa", "aligned.aln",
                        find_executable(args.m))
    alignment = mut.AlignIO.read("blast_record.aln", "clustal")
    # compute zMIc matrix
    edited = mut.prune_id_gaps(alignment, "target")
    gapped_list = mut.get_all_gaps(edited)
    edited = mut.prune(edited, gapped_list)
    (minlist, maxlist) = mut.get_extreme_columns(edited, args.low, args.high,
                                                 20)
    edited = mut.prune(edited, minlist + maxlist)
    # compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
    MI_matrix = mut.mutual_info_matrix(edited, args.b)
    ncps_array = mut.NCPS_matrix(edited, args.b)
    MIc_matrix = MI_matrix - ncps_array
    zMIc_matrix = mut.standardise_matrix(MIc_matrix)
    # check hits produced and precision compared to
    cutoff_list = []
    hit_list = []
    precision_list = []
    for cutoff in np.linspace(0.0, 3.0, num=60):
        cutoff_list.append(cutoff)
        tmatrix = mut.get_level_matrix(zMIc_matrix, cutoff)
        hits = matrix_hits(tmatrix)
        hit_list.append(hits)
        cm_residue_pairs = mut.retrieve_residue_positions(tmatrix, gapped_list,
                                                          minlist + maxlist)
        print(cm_residue_pairs)
        count = 0
        for rp in cm_residue_pairs:
            try:
                count += cont_matrix[rp[0], rp[1]]
            except IndexError:
                pass
        precision_list.append(count/hits)
    # produce the plot
    plot_twin_curves(cutoff_list, hit_list, precision_list)
