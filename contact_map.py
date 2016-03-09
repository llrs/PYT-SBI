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
import matplotlib.cm as cm
from Bio.PDB.PDBParser import PDBParser


def calc_residue_dist(residue_one, residue_two, atom):
    """Returns the distance between two residues between the selected atom."""
    logging.debug("Calculating distance between {} and {} at {} atom.".format(
                 residue_one.id, residue_two.id, atom))
    distance = residue_one[atom] - residue_two[atom]
    return distance


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


def main(structure, atom):
    """Creates the distance map between standard amino acids of a given pdb."""
    logging.debug("Reading residues of structure {}".format(structure))
    residues = tuple(structure.get_residues())
    # Filter those who are not an aminoacid

    logging.info("Filtering non-standard residues")
    residues = tuple(filter(lambda x: x.id[0] == " ", residues))
    logging.debug("Remaining {} residues.".format(len(residues)))
    dist_matrix = calc_dist_matrix(residues, atom)
    # TODO: Do a properly system log
    logging.info("Minimum distance {}".format(np.min(dist_matrix)))
    logging.info("Maximum distance {}".format(np.max(dist_matrix)))
    return dist_matrix


def contact_map(distance_map, atom):
    """Given a distance map select which residues have a relevant contact."""
    logging.debug("Selecting the contacts between atoms based on {}.".format(
                                                                    atom))
    sizes = {"CA": 15, "CB": 12}
    size = len(distance_map)
    answer = np.zeros((size, size), str)
    contact = 0
    for c in range(size):
        for b in range(size):
            # To avoid marking as contacting atoms those who simply are
            # Close in the sequence but to catch Beta-turns which are
            # of 4 residues the minimum distance is 3
            if abs(c-b) <= 2:
                answer[c][b] = " "
            elif dist_map[b][c] <= sizes[atom]:
                contact += 1
                answer[c][b] = "*"
            else:
                answer[c][b] = " "
    logging.info("Found {} contacts between residues.".format(contact/2))
    return(answer)

if __name__ == "__main__":
    logging.basicConfig(filename='contact_map.log',level=logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    msg = 'A module that calculates distance map'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a",
                           help="""Atom to calculate distance with\n
                           CA: Carbon Alpha, CB: Carbon Beta""",
                           default="CA", choices=["CA", "CB"])

    args = argparser.parse_args()
    parser = PDBParser(PERMISSIVE=1)
    logging.captureWarnings(True)
    structure = parser.get_structure("test", args.file)

    dist_map = main(structure, args.a)
    contact_map(dist_map, args.a)
    plt.imshow(dist_map, interpolation='none')
    heatmap = plt.pcolormesh(dist_map)
#     heatmap.ax.set_xlim(len(dist_map))
#     heatmap.ax.set_ylim(len(dist_map))
    plt.title('Heat map of the file {}'.format(args.file))
    legend = plt.colorbar(heatmap)
    legend.set_label("Angstroms")
    plt.savefig('distance_map_{}.png'.format(args.file))
    logging.captureWarnings(False)
