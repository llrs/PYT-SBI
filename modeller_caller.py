#!/usr/bin/env python
# encoding: utf-8

'''
Provides functions to modellize the structure of a protein.
Created on Mar 12, 2016

@author: Lluís
'''

import argparse
import logging

import matplotlib.pyplot as plt
import modeller.scripts
import modeller.automodel
from Bio.PDB import PDBList


def pdb_download(code, path=None):
    """Downloads the structure of the pdb on a file.

    Returns the file name where it is stored"""
    pdbl = PDBList()
    if path is None:
        file = pdbl.retrieve_pdb_file('1FAT')
    else:
        file = pdbl.retrieve_pdb_file('1FAT', pdir=path)
    return file

env = modeller.environ()


class modeller_caller(object):
    """Class that gets the environment variables"""
    def __init__(self, env):
        """Set the environment variables as a property"""
        self.env = env

    def convert_ali(self, fasta, pir):
        """An alignment in fasta format is converted to Modeller/pir format."""
        file_ali = "file_alignment.fasta"
        aln = modeller.alignment(self.env)
        aln.append(file=file_ali, alignment_format="FASTA", remove_gaps=False)
        aln.write(file=pir, alignment_format='PIR')
        # aln.check()  # Not sure if needed see the link:
        # https://salilab.org/modeller/9v2/manual/node269.html

    def asses_energy(self, pdb_file, name=None):
        """Asses energy of a pdb.

        Returns a matrix of energy for a plot.
        pdb_file is the name of the file to analyze
        name is the name of the output file with the energy profile"""
        logging.captureWarnings(True)
        env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
        env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters

        mdl = modeller.scripts.complete_pdb(self.env, pdb_file,
                            model_segment=('FIRST:A', 'LAST:Z'))
        s = modeller.selection(mdl)
        output_op = 'ENERGY_PROFILE NO_REPORT VERY_LONG GRADIENT'
        if name is not None:
            score = s.assess_normalized_dope()
#             score = s.assess_dope(output=output_op,
#                         file='{}.profile'.format(name),
#                         normalize_profile=True, smoothing_window=15)
        else:
            score = s.assess_normalized_dope()
#             score = s.assess_dope(output=output_op,
#                           normalize_profile=True, smoothing_window=15)
        logging.captureWarnings(False)
        return score

    # Read it and create a model
    # a = modeller.automodel.automodel(env, alnfile='output.pir',
    #         knowns='1bdmA', sequence='TvLDH',
    #         assess_methods=(modeller.automodel.assess.DOPE,
    #         #soap_protein_od.Scorer(),
    #         modeller.automodel.assess.GA341))
    #
    # a.starting_model = 1
    # a.ending_model = 5
    # a.make()

modeler = modeller_caller(env)
energies = modeler.asses_energy("pdb1cd8.ent", "profile_1cd8")
print(energies)


def plot_energy(energy):
    """Creates a plot of the energies of the pdb structure.

    Several energies can be provided and will be superposed on the plot."""
#     print("doing something with energies")
    # Plot energies of the model
    plt.plot(energy, linewidth=2)
    plt.show()