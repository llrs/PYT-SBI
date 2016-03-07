#!/usr/bin/env python
# encoding: utf-8

'''
Created on Mar 7, 2016

@author: Lluís, Leo, Ferran
'''
import argparse
import contact_map as cm


if __name__ == '__main__':
    msg = 'Predicts correlation between mutation and sturcture'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a", "--atom",
                           help="Atom to calculate distance with.",
                           default="CA", choices=["CA", "CB", "N"])
    args = argparser.parse_args()

    cm.main(args.file, args.atom)
