#!/usr/bin/env python
# encoding: utf-8

'''
Created on Mar 7, 2016

@author: Lluís, Leo, Ferran
From: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/
     python/protein_contact_map/
'''
import numpy as np
import argparse
import logging
import matplotlib.pyplot as plt
import os

from Bio.PDB.PDBParser import PDBParser


def calc_residue_dist(residue_one, residue_two, atom):
    """Returns the distance between two residues between the selected atom."""
    logging.debug("Calculating distance between {} and {} at {} atom.".format(
                 residue_one.id, residue_two.id, atom))
    distance = residue_one[atom] - residue_two[atom]
    return distance


def calc_min_dist(residue_one, residue_two):
    """Returns the minimum distance between two residues."""
    logging.debug("Calculating minimum distance between {} and {}.".format(
                 residue_one.id, residue_two.id))
    distances = []
    for a1 in residue_one:
        for a2 in residue_two:
            distances.append(a1-a2)

    return min(distances)


def comp_dist(residue_one, residue_two):
    """Compares if there is any difference distances"""
    logging.debug("comparingdistance between {} and {}.".format(
                 residue_one.id, residue_two.id))
    CA_dist = calc_residue_dist(residue_one, residue_two, "CA")
    min_dist = calc_min_dist(residue_one, residue_two)
    if CA_dist > 12 and min_dist < 6:
        return 1
    else:
        return 0


def calc_dist_matrix(chain, atom):
    """Returns a matrix of distances between two chains."""
    logging.debug("Calculating distance matrix for {}".format(chain))
    size = len(chain)
    answer = np.zeros((size, size), np.float)
    for row, residue_one in enumerate(chain):
        for col, residue_two in enumerate(chain):
            answer[row, col] = calc_residue_dist(residue_one, residue_two,
                                                 atom)
    return answer


def calc_mdist_matrix(chain):
    """Returns a matrix of distances between two chains."""
    logging.debug("Calculating distance matrix for {}".format(chain))
    size = len(chain)
    answer = np.zeros((size, size), np.float)
    for row, residue_one in enumerate(chain):
        for col, residue_two in enumerate(chain):
            answer[row, col] = calc_min_dist(residue_one, residue_two)
    return answer


def calc_matrix(chain):
    """Calculate those who happen to meet the threeshold """
    logging.debug("Calculating distance matrix for {}".format(chain))
    size = len(chain)
    answer = np.zeros((size, size), np.float)
    for row, residue_one in enumerate(chain):
        for col, residue_two in enumerate(chain):
            distance = []
            for atom1 in residue_one:
                for atom2 in residue_two:
                    distance.append(atom1 - atom2)
            answer[row, col] = calc_min_dist(residue_one, residue_two)
    return answer


def calc_relevant(chain):
    """Count how many non-identified positions by CA"""
    val = sum([comp_dist(residue_one, residue_two) for residue_one in chain
               for residue_two in chain])
    print(val)


def main(structure, atom=None):
    """Creates the distance map between standard amino acids of a given pdb."""
    logging.debug("Reading residues of structure {}".format(structure))
    residues = tuple(structure.get_residues())
    # Filter those who are not an aminoacid

    logging.info("Filtering non-standard residues")
    residues = tuple(filter(lambda x: x.id[0] == " ", residues))
    logging.debug("Remaining {} residues.".format(len(residues)))

#     if atom is None:
#         dist_matrix = calc_mdist_matrix(residues)
#     else:
#         dist_matrix = calc_dist_matrix(residues, atom)
    calc_relevant(residues)
    # TODO: Do a properly system log
#     logging.info("Minimum distance {}".format(np.min(dist_matrix)))
#     logging.info("Maximum distance {}".format(np.max(dist_matrix)))
#     return dist_matrix


def contact_map(distance_map, atom):
    """Given a distance map select which residues have a relevant contact."""
    logging.debug("Selecting the contacts between atoms based on {}.".format(
                                                                    atom))
    sizes = {"CA": 15, "CB": 12, None: 6}
    size = len(distance_map)
    answer = [[False for j in range(size)] for i in range(size)]
    contact = 0
    for c in range(size):
        for b in range(size):
            # To avoid marking as contacting atoms those who simply are
            # Close in the sequence but to catch Beta-turns which are
            # of 4 residues the minimum distance is 3
            if abs(c-b) <= 2:
                pass
            elif dist_map[c][b] <= sizes[atom]:
                contact += 1
                answer[c][b] = True #dist_map[c][b]
#             else:
#                 answer[c][b] = 30
    logging.info("Found {} contacts between residues.".format(contact/2))
    return(answer)

if __name__ == "__main__":

    logging.basicConfig(filename='contact_map.log', level=logging.DEBUG)
    fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
     - %(message)s"""
    formatter = logging.Formatter(fmt)
    msg = 'A module that calculates distance map'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a",
                           help="""Atom to calculate distance with
                           CA: Carbon Alpha, CB: Carbon Beta. By default is
                            the minimum distance between residues""",
                           default=None, choices=["CA", "CB"])

    args = argparser.parse_args()
    base = os.path.basename(args.file)
    name_f = os.path.splitext(base)[0]
    parser = PDBParser(PERMISSIVE=1)
    logging.captureWarnings(True)
    structure = parser.get_structure("test", args.file)

    dist_map = main(structure, args.a)
#     cont_map = contact_map(dist_map, args.a)
# 
#     # Plot the distance map
#     plt.imshow(dist_map, interpolation='none')
#     heatmap = plt.pcolormesh(dist_map)
#     plt.title('Heat map of the file {}'.format(name_f))
#     legend = plt.colorbar(heatmap)
#     legend.set_label("Angstroms")
#     plt.savefig('distance_map_{}_{}.png'.format(name_f, args.a))
# 
#     logging.captureWarnings(False)
# 
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
# 
#     ax.imshow(cont_map, aspect='auto',
#               cmap=plt.cm.gray, interpolation='nearest')
#     fig.savefig("contact_map_{}_{}.png".format(name_f, args.a))
#     fig.show()
