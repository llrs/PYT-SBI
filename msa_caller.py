#!/usr/bin/python3
"""
Provides a means to call 3 different MSA methods
Created on Mar 13, 2016
@author: Leo, Ferran, Lluís
"""
# standard modules
import os
import argparse
from distutils.spawn import find_executable

# Biopython modules
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import TCoffeeCommandline

def call_msa_method(method_label, in_file, out_file, executable_path):
	"""Calls the appropriate MSA method"""
	if method_label == "clustalw":
		cline = ClustalwCommandline(executable_path, 
								    infile = in_file)
	elif method_label == "muscle":
		cline = MuscleCommandline(executable_path, 
								  input = in_file, 
								  out = out_file)
	elif method_label == "t_coffee":
		cline = TCoffeeCommandline(infile = in_file,
								   output = "clustalw",
								   outfile = out_file)
	else:
		return
	assert os.path.isfile(executable_path), "Executable missing"
	stdout, stderr = cline()

if __name__ == '__main__':
	msg = 'Carries out a MSA'
	argparser = argparse.ArgumentParser(description=msg)
	argparser.add_argument("-f",
						   help = "File with sequences to be aligned.")
	argparser.add_argument("-m", 
                           help = "Method for MSA.",
                           default = "clustalw")
	argparser.add_argument("-o", help = "Output file.")
	args = argparser.parse_args()
	call_msa_method(args.m, args.f, "aligned.aln", find_executable(args.m))
