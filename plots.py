#!/usr/bin/python3
"""
Deals with plots and graphical outputs.
Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import matplotlib.pyplot as plt
import numpy as np
import logging

# non-standard modules
import mutual_information as mut


def plot_heatmap(distances, name_file, title, option):
    """Plots the distances between the residues."""
    logging.info("Plotting the distance map for {}".format(name_file))
    plt.imshow(distances, interpolation='none')
    heatmap = plt.pcolormesh(distances)
    plt.title(title)
    legend = plt.colorbar(heatmap)
    legend.set_label("Angstroms")
    plt.savefig('heatmap_{}_{}.png'.format(name_file, option),
                format="png")

def plot_matrix_binary(matrix, name_file, title, option):
    """Plots the contact map between residues black means contact."""
    logging.info("Plotting the contact map for {}".format(name_file))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(title)
    imgplot = ax.imshow(matrix, cmap='Greys', interpolation='none')
    plt.savefig('contact_map_{}_{}.png'.format(name_file, option),
                format="png")


def plot_twin_curves(cutoff_list, hit_list, precision_list, name_file):
    """Plots the precision and the hits depending on the cutoff."""
    fig, ax1 = plt.subplots()
    ax1.plot(cutoff_list, hit_list, 'ob-')
    ax1.set_xlabel('Cutoff level L')
    # Make the y-axis label match the line color.
    ax1.set_ylabel('Predicted CM residue pairs', color='b')
    # Make the second plot with another axis
    ax2 = ax1.twinx()
    ax2.plot(cutoff_list, precision_list, 'or-')
    ax2.set_ylabel('Precision', color='r')
    plt.savefig('cutoff_zMIc_{}.png'.format(name_file),
                format="png")
#     plt.show()


def precision_analysis(zMIc_matrix, cont_matrix, gapped_list, mm, l=0.0, h=3.0, num=60):
    """Calculates precision and hits for different threshold levels of zMIc.

   zMIc matrix: with the score of zMIc
   cont_matrix: matrix of contact residues by a default distance.
   gapped_list: columns of the zMIc that had gaps and weren't used
   mm list the sum of minlist and maxlist with prunned columns
   l, h, num: are the parameters for np.linspace
   l: start, h:stop, num: the number of intervals to generate

   Return a list with the cutoff, the number of hits, and the precision foreach
   interval"""
    cutoff_list = []
    hit_list = []
    precision_list = []
    for cutoff in np.linspace(l, h, num):
        cutoff_list.append(cutoff)
        tmatrix = mut.get_level_matrix(zMIc_matrix, cutoff)
        hits = mut.matrix_hits(tmatrix)
        hit_list.append(hits)
        cm_residue_pairs = mut.retrieve_residue_positions(tmatrix, gapped_list,
                                                          mm)
        print(cm_residue_pairs)
        count = 0
        for rp in cm_residue_pairs:
            count += cont_matrix[rp[0], rp[1]]
        precision_list.append(count/hits)
    else:
        return(cutoff_list, hit_list, precision_list)


if __name__ == "__main__":
    pass
