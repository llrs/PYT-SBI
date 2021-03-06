#!/usr/bin/python3
# encoding: utf-8
"""
Provides a means to call 3 different MSA methods
Created on Mar 13, 2016
@author: Leo, Ferran, Lluís
"""
# standard modules
import argparse
import logging
from distutils.spawn import find_executable

# Biopython modules
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import TCoffeeCommandline
from Bio.Application import ApplicationError


def call_msa_method(method, in_file, out_file, output_format=None):
    """Calls the appropriate program to generate a MSA.
    the muscle is the only one that only produce FASTA alignments,
    in the other programs this can be changed"""
    logging.info("Creating a MSA  with {} from {} to {}".format(
                                            method, in_file, out_file))

    if find_executable(method):
        pass
    else:
        msg = "Program {} not found it is installed?"
        raise UnboundLocalError(msg.format(method))

    if method == "clustalw":
        cline = ClustalwCommandline(method,
                                    infile=in_file, output=output_format,
                                    OUTFILE=out_file, type="PROTEIN")
    elif method == "muscle":
        cline = MuscleCommandline(method,
                                  input=in_file,
                                  out=out_file)
    elif method == "t_coffee":
        cline = TCoffeeCommandline("t_coffee",
                                   infile=in_file,
                                   output=output_format,
                                   outfile=out_file)
    try:
        stdout, stderr = cline()
    except ApplicationError:
        raise IOError("The input file doesn't contain an alignment!")
    logging.debug(stderr)
    return(stdout)

if __name__ == '__main__':
    fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
     - %(message)s"""
    logging.basicConfig(filename='msa_caller.log', level=logging.DEBUG,
                        format=fmt)

    msg = 'Carries out a MSA'
    arg_helper = argparse.ArgumentDefaultsHelpFormatter
    programs_MSA = ["clustalw", "muscle", "t_coffee"]
    argparser = argparse.ArgumentParser(description=msg,
                                        formatter_class=arg_helper)
    argparser.add_argument("file",
                           help="File with sequences to be aligned.")
    argparser.add_argument("-m",
                           help="Method for MSA.",
                           default="clustalw", choices=programs_MSA)
    argparser.add_argument("-o", help="Output file.", default="aligned.aln")
    argparser.add_argument("-f", help="Output format.", default="fasta")
    args = argparser.parse_args()
    call_msa_method(args.m, args.file, args.o, args.f)
    logging.info("Ended program")
    exit(0)
