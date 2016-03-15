#!/usr/bin/env python
# encoding: utf-8

'''
Provides functions to modellize the structure of a protein.
Created on Mar 12, 2016

@author: Lluís
'''

import argparse
import logging
import os

from Bio import AlignIO

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
        """An alignment in fasta format is converted to Modeller/pir format.

        It downloads the pdb of the alignment to fetch the necessary data for 
        the pir format."""
        assert pir != "output.pir"  # Assumption
        logging.captureWarnings(True)
        aln = modeller.alignment(self.env)
        aln.append(file=fasta, alignment_format="FASTA", remove_gaps=False)
        aln.write(file="output.pir", alignment_format='PIR')

        alignment = AlignIO.read(fasta, "fasta")
        values = []
        for record in alignment:
            values.append([len(record), record.id])
            pdb_download(record.id)  # Download the pdb to build the model
        


        self.pir = pir  # Set the pir as an attribute
        # Conver the pir into a understandable pir format?
        with open(pir, "w") as out:
            with open("output.pir", "r") as fl:
                records = fl.read()
                records = records.split(">")
                for n, record in enumerate(records, start=-1):
                    lines = record.split("\n")
                    if lines[0] == "":
                        continue
                    id_pdb = lines[0].split(";")[1]
                    lines[0] = ">"+lines[0]
                    fields = lines[1].split(":")
                    fields[1] = id_pdb
                    fields[2] = "1"
                    if values[n][1] == id_pdb.rstrip():
                        fields[3] = str(values[n][0])
                    else:
                        fields[3] = "500"  # Default length of the sequence
                    lines[1] = ":".join(fields)
                    lines_o = "\n".join(lines)
                    out.write(lines_o)
        os.remove("output.pir")
        logging.captureWarnings(False)

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

    def modelize(self, alig_pir, known, seq):
        """Uses automodel to generate a model of the """
        logging.captureWarnings(True)
        a = modeller.automodel.automodel(self.env, alnfile=alig_pir,
                knowns=known, sequence=seq,
                assess_methods=modeller.automodel.assess.DOPE)

        a.starting_model = 1
        a.ending_model = 5
        a.make()
        self.outputs = a.outputs
        logging.captureWarnings(False)

modeler = modeller_caller(env)
# modeler.convert_ali("output.fastaa", "output_modeller.pir")
modeler.modelize("output_modeller.pir", "1cd8A", "1dc8A")
# energies = modeler.asses_energy("pdb1cd8.ent", "profile_1cd8")
# print(energies)


def plot_energy(energy):
    """Creates a plot of the energies of the pdb structure.

    Several energies can be provided and will be superposed on the plot."""
#     print("doing something with energies")
    # Plot energies of the model
    plt.plot(energy, linewidth=2)
    plt.show()
