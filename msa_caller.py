#!/usr/bin/python3
"""
Provides a means to call 3 different MSA methods
Created on Mar 13, 2016
@author: Leo, Ferran, Lluís
"""
# standard modules
import os
import argparse
import logging
from distutils.spawn import find_executable

# Biopython modules
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import TCoffeeCommandline


def call_msa_method(method, in_file, out_file, executable_path, output_format):
    """Calls the appropriate MSA method"""
    logging.info("Creating a MSA  with {} from {} to {}".format(
                                            method, in_file, out_file))
    assert os.path.isfile(executable_path), "Executable missing"
    if method == "clustalw":
        cline = ClustalwCommandline(executable_path,
                                    infile=in_file, output_format="fasta")
    elif method == "muscle":
        cline = MuscleCommandline(executable_path,
                                  input=in_file,
                                  out=out_file)
    elif method == "t_coffee":
        cline = TCoffeeCommandline(infile=in_file,
                                   output="clustalw",
                                   outfile=out_file)

    stdout, stderr = cline()
    logging.debug(stderr)
    return(stdout)

if __name__ == '__main__':
    msg = 'Carries out a MSA'
    programs_MSA = ["clustalw", "muscle", "t_coffee"]
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file",
                           help="File with sequences to be aligned.")
    argparser.add_argument("-m",
                           help="Method for MSA.",
                           default="clustalw", choices=programs_MSA)
    argparser.add_argument("-o", help="Output file.")
    args = argparser.parse_args()
    call_msa_method(args.m, args.f, "aligned.aln", find_executable(args.m))
