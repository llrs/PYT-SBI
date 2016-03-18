#!/usr/bin/env python
# encoding: utf-8

'''
Created on Mar 7, 2016
Provides functions to create a contact map given a pdb structure.
@author: Lluís, Leo, Ferran
Adapted from: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock
     /python/protein_contact_map/
'''
import argparse
import logging
import os

import numpy as np
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt

import plots


def calc_residue_dist(residue_one, residue_two, atom):
    """Returns the distance between two residues between the selected atom."""
    msg = "Calculating distance between {} {} and {} {} at {} atom."
    logging.debug(msg.format(residue_one.get_resname(), residue_one.id[1],
                             residue_two.get_resname(), residue_two.id[1],
                             atom))
    assert atom is not None
    try:
        distance = residue_one[atom] - residue_two[atom]
    except KeyError:
        logging.info("Using C-Alpha distance for glycine.")
        distance = residue_one["CA"] - residue_two["CA"]
    return distance


def calc_min_dist(residue_one, residue_two):
    """Returns the minimum distance between two residues."""
    log_msg = "Calculating minimum distance between {} {} and {} {}."
    logging.debug(log_msg.format(residue_one.get_resname(), residue_one.id[1],
                                 residue_two.get_resname(), residue_two.id[1]))
    distances = []
    for a1 in residue_one:
        for a2 in residue_two:
            distances.append(a1-a2)

    return min(distances)


def calc_dist_matrix(chain, atom=None):
    """Returns a matrix of distances between residues of the same chain."""
    logging.debug("Calculating distance matrix for {}".format(chain))
    size = len(chain)
    answer = np.zeros((size, size), np.float)
    for row, residue_one in enumerate(chain):
        for col, residue_two in enumerate(chain):
            if atom == "min" or atom is None:
                answer[row, col] = calc_min_dist(residue_one, residue_two)
            else:
                answer[row, col] = calc_residue_dist(residue_one, residue_two,
                                                     atom)
    logging.info("Minimum distance {}".format(np.min(answer)))
    logging.info("Maximum distance {}".format(np.max(answer)))
    return answer


def filter_residues(structure):
    """Filters non-standard amino acids of the structure."""
    logging.debug("Reading residues of structure {}".format(structure))
    residues = tuple(structure.get_residues())
    # Filter those who are not an amino acid
    logging.debug("The structure has {} residues.".format(len(residues)))
    logging.info("Filtering non-standard residues")
    residues = tuple(filter(lambda x: x.id[0] == " ", residues))
    logging.debug("Remaining {} residues.".format(len(residues)))
    return residues


def contact_map(distance_map, atom, dist=None):
    """Given a distance map select which residues have a relevant contact."""
    logging.debug("Selecting the contacts between atoms based on {}.".format(
                                                                    atom))
    if dist is None:
        distances = {"CA": 15, "CB": 12, "min": 6}
    else:
        distances = dist

    size = len(distance_map)
    answer = np.zeros((size, size), dtype=bool)
    contact = 0
    for c in range(size):
        for b in range(size):
            # To avoid marking as contacting atoms those who simply are
            # Close in the sequence but to catch Beta-turns which are
            # of 4 residues the minimum distance is 3
            if abs(c-b) <= 2:
                pass
            elif distance_map[c][b] <= distances[atom]:
                contact += 1
                answer[c][b] = True
    logging.info("Found {} contacts between residues.".format(contact/2))
    return(answer)


def plot_distance(distances, name_file, title, option):
    """Plots the distances between the residues."""
    logging.info("Plotting the distance map for {}".format(name_file))
    plt.imshow(distances, interpolation='none')
    heatmap = plt.pcolormesh(distances)
    plt.title(title)
    legend = plt.colorbar(heatmap)
    legend.set_label("Angstroms")
    plt.savefig('distance_map_{}_{}.png'.format(name_file, option),
                format="png")


def plot_contacts(contacts, name_file, title, option):
    """Plots the contact map between residues"""
    logging.info("Plotting the contact map for {}".format(name_file))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(title)
    ax.imshow(contacts, aspect='auto',
              cmap=plt.cm.gray, interpolation='nearest')
    plt.savefig("contact_map_{}_{}.png".format(name_file, option),
                format="png")


if __name__ == "__main__":

    logging.basicConfig(filename='contact_map.log', level=logging.DEBUG)
    fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
     - %(message)s"""
    formatter = logging.Formatter(fmt)
    msg = 'A module that calculates distance map'
    argparser = argparse.ArgumentParser(description=msg,
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a",
                           help="""Atom to calculate distance with
                           CA: Carbon Alpha, CB: Carbon Beta.""",
                           default="min", choices=["CA", "CB", "min"])
    argparser.add_argument("-CA",
                           help="""Set the threshold distance for Carbon alpha
                           """,
                           type=int,
                           default=15)
    argparser.add_argument("-CB",
                           help="""Set the threshold distance for Carbon beta
                           """,
                           type=int,
                           default=12)
    argparser.add_argument("-min",
                           help="""Set the minimal threshold distance""",
                           type=int,
                           default=16)
    args = argparser.parse_args()
    dist = {"CA": args.CA, "CB": args.CB, "min": args.min}
    base = os.path.basename(args.file)
    name_f = os.path.splitext(base)[0]
    parser = PDBParser(PERMISSIVE=1)
    logging.captureWarnings(True)

    structure = parser.get_structure("test", args.file)

    residues = filter_residues(structure)
    dist_matrix = calc_dist_matrix(residues, args.a)
    title_dist = 'Distances of the file {}'.format(name_f)
    plots.plot_distance(dist_matrix, name_f, title_dist, args.a)
    cont_matrix = contact_map(dist_matrix, args.a, dist)
    title_cont = 'Contacts of the file {}'.format(name_f)
    plots.plot_contacts(cont_matrix, name_f, title_cont, args.a)
    logging.captureWarnings(False)
