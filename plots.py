#!/usr/bin/python3
# encoding: utf-8
"""
Deals with plots and graphical outputs, some calculus
Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import logging
import os

import matplotlib.pyplot as plt
import numpy as np

# non-standard modules
from Bio.PDB import PDBList

import mutual_information as mut


def plot_heatmap(matrix, name_file, title, option, label):
    """Plots a heatmap of a matrix."""

    logging.info("Plotting the distance map for {}".format(name_file))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(title)
    heatmap = plt.pcolormesh(matrix)
    legend = plt.colorbar(heatmap)
    legend.set_label(label)
    imgplot = ax.imshow(matrix, interpolation='none')
    name_out = 'heatmap_{}_{}.png'.format(name_file, option)
    plt.savefig(name_out, format="png")
    return name_out


def plot_matrix_binary(matrix, name_file, title, option):
    """Plots the contact map between residues black means contact."""

    logging.info("Plotting the contact map for {}".format(name_file))

    fig_b = plt.figure()
    ax = fig_b.add_subplot(111)
    fig_b.suptitle(title)
    imgplot = ax.imshow(matrix, cmap='Greys', interpolation='none')
    name_out = 'contact_map_{}_{}.png'.format(name_file, option)
    plt.savefig(name_out, format="png")
    return name_out


def plot_twin_curves(cutoff_list, hit_list, precision_list, name_file):
    """Plots the precision and the hits depending on the cutoff."""

    logging.info("Plotting the precision and hits against zMIc cutoff.")
    fig, ax1 = plt.subplots()
    ax1.plot(cutoff_list, hit_list, 'ob-')
    ax1.set_xlabel('Cutoff level L')
    # Make the y-axis label match the line color.
    ax1.set_ylabel('Predicted CM residue pairs', color='b')
    # Make the second plot with another axis
    ax2 = ax1.twinx()
    ax2.plot(cutoff_list, precision_list, 'or-')
    ax2.set_ylabel('Precision', color='r')
    name_out = 'cutoff_zMIc_{}.png'.format(name_file)
    plt.savefig(name_out, format="png")
    return name_out
#     plt.show()


def precision_analysis(zMIc_m, cont_m, gapped_list, mm, l=0.0, h=3.0, num=60):
    """Calculates precision and hits for different threshold levels of zMIc.

    zMIc_m: matrix with the score of zMIc
    cont_m: matrix of contact residues by a default distance.
    gapped_list: columns of the zMIc that had gaps and weren't used
    mm list the sum of minlist and maxlist with prunned columns
    l, h, num: are the parameters for np.linspace
    l: start, h:stop, num: the number of intervals to generate

    Return a list with the cutoff, the number of hits and the precision for
    each interval"""

    logging.info("Calculating the precision and number of hits.")

    cutoff_list = []
    hit_list = []
    precision_list = []
    for cutoff in np.linspace(l, h, num):
        cutoff_list.append(cutoff)
        tmatrix = mut.get_level_matrix(zMIc_m, cutoff)
        hits = mut.matrix_hits(tmatrix)
        hit_list.append(hits)
        cm_residue_pairs = mut.retrieve_residue_positions(tmatrix, gapped_list,
                                                          mm)
#         print(cm_residue_pairs)
        count = 0
        for rp in cm_residue_pairs:
            count += cont_m[rp[0], rp[1]]
        precision_list.append(count/hits)
    else:
        return(cutoff_list, hit_list, precision_list)


def pdb_download(code, path=None):
    """Downloads the structure of the pdb on a file.

    cod is the pdb code of the structure
    path is the localization where it will be downloaded

    Returns the file name where it is stored"""

    logging.info("Downloading pdb %s.", code)

    logging.captureWarnings(True)
    pdbl = PDBList(obsolete_pdb=os.getcwd())
    if path is None:
        file = pdbl.retrieve_pdb_file(code)
    else:
        file = pdbl.retrieve_pdb_file(code, pdir=path)
    logging.captureWarnings(False)
    return file


def distances_zMIcs(dict_pairs, dist_matrix):
    """Relates the distance_matrix with the zMIc matrix."""

    logging.info("Relating distances with zMIc values.")
    print(dict_pairs)
    list_coor = dict_pairs.keys()
    list_zMIc = dict_pairs.values()
    list_dist = []
#     for i in range(len(list_coor)):
#         list_dist.append(dist_matrix[list_coor[i][0]][list_coor[i][1]])
    for element in list_coor:
        list_dist.append(dist_matrix[element[0]][element[1]])
    assert len(list_dist) == len(list_zMIc)
    return (list_dist, list_zMIc)


def plot_distance_zMIc(list_dist, list_zMIc, name_file):
    """Plots the zMIc values of residues pairs against the distance."""

    msg = """Plotting the zMIc values of residues pairs facing the distance
    between them."""
    logging.info(msg)

    fig = plt.figure()
    fig_title = 'zMIc  VS. distance between residues pairs of {}'
    fig.suptitle(fig_title.format(name_file))
    assert len(list_dist) == len(list_zMIc)
    plt.scatter(list_zMIc, list_dist, color='blue', s=5, edgecolor='none')
    plt.xlabel('zMIc values (z-scores)', fontsize=16)
    plt.ylabel('Distances (Angstroms)', fontsize=16)
    file_out = 'distance_zMIc_{}.png'.format(name_file)
    plt.savefig(file_out)
    return file_out

if __name__ == "__main__":
    print("""This module provides some functions used by several part of the
    program.\nIt should not be run as an independent module.""")
