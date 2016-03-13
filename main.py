#!/usr/bin/env python
# encoding: utf-8

'''
Created on Mar 7, 2016

@author: Lluís, Leo, Ferran
'''
import argparse
import contact_map as cm


## compare the maps using CA and minimum distance to search for something.


if __name__ == '__main__':
    msg = 'Predicts correlation between mutation and structure'
    argparser = argparse.ArgumentParser(description=msg)
    argparser.add_argument("file", help="PDB structure to analyze.")
    argparser.add_argument("-a",
                           help="Atom to calculate distance with.",
                           default=None, choices=["CA", "CB"])
    args = argparser.parse_args()
