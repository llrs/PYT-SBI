#!/usr/bin/env python
# encoding: utf-8

'''
Created on Mar 7, 2016

@author: Lluís, Leo, Ferran
From: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/
     python/protein_contact_map/
'''

from Bio.PDB.PDBParser import PDBParser
import numpy as np
import argparse


def calc_residue_dist(residue_one, residue_two, atom):
    """Returns the distance between two residues between the selected atom."""
    distance = residue_one[atom] - residue_two[atom]
    return distance


def calc_dist_matrix(chain, atom):
    """Returns a matrix of distances between two chains."""
    size = len(chain)
    answer = np.zeros((size, size), np.float)
    for row, residue_one in enumerate(chain):
        for col, residue_two in enumerate(chain):
            answer[row, col] = calc_residue_dist(residue_one, residue_two,
                                                 atom)
    return answer


def main(file, atom):
    """Creates the distance map between aminoacids of a given pdb."""
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("test", file)
    # for model in structure:

    residues = tuple(structure.get_residues())

    # Filter those who are not an aminoacid
    residues = tuple(filter(lambda x: x.id[0] == " ", residues))

    dist_matrix = calc_dist_matrix(residues, atom)
    # TODO: Do a properly system log
    print("Minimum distance", np.min(dist_matrix))
    print("Maximum distance", np.max(dist_matrix))
    return dist_matrix

if __name__ == "__main__":
    msg = 'A module that calculates distance map'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a", "--atom",
                           help="Atom to calculate distance with.",
                           default="CA", choices=["CA", "CB", "N"])
    args = argparser.parse_args()
    main(args.file, args.atom)
