#!/usr/bin/env python
# encoding: utf-8

'''
Provides functions to modellize the structure of a protein. 
Created on Mar 12, 2016

@author: Lluís
'''

import argparse

import matplotlib as plt
import modeller.scripts
import modeller.automodel
from Bio.PDB import PDBList


env = modeller.environ()
# Read a fasta alingment and prepare for a model
# file_ali = "file_alignment.fasta"
# aln = modeller.alignment(env)
# aln.append(file=file_ali, alignment_format="FASTA", remove_gaps=False)
# aln.write(file='output.pir', alignment_format='PIR')
# aln.check()  # Not sure if needed see the link:
# https://salilab.org/modeller/9v2/manual/node269.html
# Download the pdb file
# pdbl = PDBList()
# file = pdbl.retrieve_pdb_file('1FAT', pdir="path to the download folder")
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


# Asses energy profile
def asses_energy(pdb_file):
    """"""
    env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters

    mdl = modeller.scripts.complete_pdb(env, 'pdb1cd8.ent',
                        model_segment=('FIRST:A', 'LAST:Z'))  # Model of pdb
    s = modeller.selection(mdl)
    a = s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='1cd8.profile',
              normalize_profile=True, smoothing_window=15)
    return a
# Plot energies of the model
# plt.figure(energies)
# plt.show()
