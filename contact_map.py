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

import plots


fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
     - %(message)s"""

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


def main(file, atom, CAd=15, CBd=12, mind=6):
    """Analyze the pdb using distance between atom and minimum distances."""

    logging.info("Analyzing %s using %s", file, atom)

    dist = {"CA": CAd, "CB": CBd, "min": mind}
    base = os.path.basename(args.file)
    name_f = os.path.splitext(base)[0]
    parser = PDBParser(PERMISSIVE=1)
    logging.captureWarnings(True)

    structure = parser.get_structure("test", file)

    residues = filter_residues(structure)
    dist_matrix = calc_dist_matrix(residues, atom)
    title_dist = 'Distances of the file {}'.format(name_f)
    name_heatmap = plots.plot_heatmap(dist_matrix, name_f, title_dist, atom)
    logging.info("Heatmap %s created", name_heatmap)
    cont_matrix = contact_map(dist_matrix, atom, dist)
    title_cont = 'Contacts of the file {}'.format(name_f)
    name_bin = plots.plot_matrix_binary(cont_matrix, name_f, title_cont, atom)
    logging.info("Contact map %s created", name_bin)
    logging.captureWarnings(False)

    return(dist_matrix, cont_matrix)

if __name__ == "__main__":

    logging.basicConfig(filename='contact_map.log', level=logging.DEBUG)
    
    formatter = logging.Formatter(fmt)
    msg = 'A module that calculates distance map'
    args_helper = argparse.ArgumentDefaultsHelpFormatter
    argparser = argparse.ArgumentParser(description=msg,
                                        formatter_class=args_helper)
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
                           default=6)
    args = argparser.parse_args()

    main(args.file, args.a, args.CA, args.CB, args.min)
#     dist = {"CA": args.CA, "CB": args.CB, "min": args.min}
#     base = os.path.basename(args.file)
#     name_f = os.path.splitext(base)[0]
#     parser = PDBParser(PERMISSIVE=1)
#     logging.captureWarnings(True)
# 
#     structure = parser.get_structure("test", args.file)
# 
#     residues = filter_residues(structure)
#     dist_matrix = calc_dist_matrix(residues, args.a)
#     title_dist = 'Distances of the file {}'.format(name_f)
#     plots.plot_heatmap(dist_matrix, name_f, title_dist, args.a)
#     cont_matrix = contact_map(dist_matrix, args.a, dist)
#     title_cont = 'Contacts of the file {}'.format(name_f)
#     plots.plot_matrix_binary(cont_matrix, name_f, title_cont, args.a)
#     logging.captureWarnings(False)
