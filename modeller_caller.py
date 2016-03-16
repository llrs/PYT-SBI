#!/usr/bin/python3
# encoding: utf-8

'''
Provides functions to modellize the structure of a protein.
Created on Mar 12, 2016

@author: Lluís
'''

import argparse
import logging
import os
import ftplib
import urllib

from Bio import SeqIO

import matplotlib.pyplot as plt
import modeller.scripts
import modeller.automodel
from Bio.PDB import PDBList

env = modeller.environ()  # Some variables needed for the modeller


def pdb_download(code, path=None):
    """Downloads the structure of the pdb on a file.

    Returns the file name where it is stored"""
    logging.info("Downloading pdb %s.", code)
    pdbl = PDBList(obsolete_pdb=os.getcwd())
    if path is None:
        file = pdbl.retrieve_pdb_file(code)
    else:
        file = pdbl.retrieve_pdb_file(code, pdir=path)
    return file


class modeller_caller(object):
    """Class that gets the environment variables"""

    def __init__(self, env):
        """Set the environment variables as a property of the class."""
        self.env = env

    def _extract_id(self, header):
        """Returns the id of the sequence."""
        return header.split("|")[-2]

    def convert_ali(self, fasta, pir):
        """An alignment in fasta format is converted to Modeller/pir format.

        It downloads the pdb of the alignment to fetch the necessary data for
        the pir format."""
        assert pir != "output.pir"  # Assumption
        logging.captureWarnings(True)
        aln = modeller.alignment(self.env)
        aln.append(file=fasta, alignment_format="FASTA", remove_gaps=False)
        aln.write(file="output.pir", alignment_format='PIR')
        fasta_h = open(fasta, "r")
        sequences = SeqIO.parse(fasta_h, "fasta")
        values = []
        for record in sequences:
            values.append([len(record), self._extract_id(record.id)])
#             print(record.id)
            # Download the pdb to build the model
            # modeller search for all the posible names of the file
            try:
                pdb_download(self._extract_id(record.id), os.getcwd())
            except urllib.error.URLError:
                pass
            except ftplib.error_perm:
                pass

        self.pir = pir  # Set the pir as an attribute
        # Convert the pir into a understandable pir format?
        with open(pir, "w") as out:
            with open("output.pir", "r") as fl:
                records = fl.read()
                records = records.split(">")
                for n, record in enumerate(records, start=-1):
                    lines = record.split("\n")
                    if lines[0] == "":
                        continue
                    id_pdb = self._extract_id(lines[0])
                    lines[0] = ">"+lines[0].split(";")[0]+";"+id_pdb
                    fields = lines[1].split(":")
                    fields[0] = "structureX"
                    fields[1] = id_pdb
                    fields[2] = "1"
                    if values[n][1] == id_pdb.rstrip():
                        fields[4] = str(values[n][0])
                    else:
                        fields[4] = "500"  # Default length of the sequence
                    fields_a = []
                    for field in fields[:]:
                        if field == "":
                            fields_a.append(".")
                        else:
                            fields_a.append(field)
                    lines[1] = ":".join(fields_a)
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
        # Setting to read all the models (not sure if it is the right way)
        mdl = modeller.scripts.complete_pdb(self.env, pdb_file,
                            model_segment=('FIRST:A', 'LAST:Z'))
        s = modeller.selection(mdl)
        output_op = 'ENERGY_PROFILE NO_REPORT VERY_LONG GRADIENT'
        if name is not None:
            score = s.assess_dope(output=output_op,
                        file='{}.profile'.format(name),
                        normalize_profile=True, smoothing_window=15)
        else:
            score = s.assess_normalized_dope(normalize_profile=True,
                                             smoothing_window=15)
        logging.captureWarnings(False)
        return score

    def modelize(self, alig_pir, known, seq):
        """Uses automodel to generate a model of the protein.

        alig_pir is the alignment in pir format of the proteins.
        known are the pdb structures known which are similar to the sequence.
        it can be a list of pdb id which should be on the same folder and
        in the alignment file
        seq is the sequence we want to create the structure"""
        logging.captureWarnings(True)
        a = modeller.automodel.automodel(self.env, alnfile=alig_pir,
                knowns=known, sequence=seq,
                assess_methods=modeller.automodel.assess.DOPE)

        a.starting_model = 1
        a.ending_model = 5
        a.make()
        self.outputs = a.outputs
        logging.captureWarnings(False)
        return self.outputs


def plot_energy(energy):
    """Creates a plot of the energies of the pdb structure.

    Several energies can be provided and will be superposed on the plot."""
#     print("doing something with energies")
    # Plot energies of the model
    plt.plot(energy, linewidth=2)
    plt.show()

if __name__ == "__main__":
    msg = 'Creates models of the sequences'
    argparser = argparse.ArgumentParser(description=msg,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("seq", help="Name of the sequence to analyse")
    argparser.add_argument("models", help="Models of the file")
    argparser.add_argument("-pir", help="Name of the file in pir format")
    argparser.add_argument("-fasta",
                           help="File with sequences in fasta format")
    args = argparser.parse_args()
#     print(os.getcwd())
#     pdb_download("3kud", os.getcwd())
    env = modeller.environ()  # Some variables needed for the modeller
    modeler = modeller_caller(env)

    # Convert the fasta alignment in pir format
    if not args.fasta and not args.pir:
        raise argparser.error("Required a fasta or a pir alignment")
    elif args.fasta:
        modeler.convert_ali(args.fasta, args.pir)
    # modeler.convert_ali("output.fastaa", "output_modeller.pir")
    else:
        modeler.modelize(args.pir, args.seq, args.models)
