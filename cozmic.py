#!/usr/bin/pyhton3
# encoding: utf-8
'''
Deals with the main workflow of the program. Given a protein, provides
contact information either from PDB structural data or by comparative 
modelling. It also provides full computations to analyse correlated 
mutation signals, including search for homologs, MSA and computation of 
MIc Z-scores. Finally, it plots graphics for visualisation and analysis 
of the associations between contact, residue distances and CM signal 
in the form of MIc Z-scores.
Created on Mar 7, 2016
@author: Leo Madsen, Ferran Muiños, Lluís Revilla Sancho
'''
# Standard modules
import argparse
import matplotlib.pyplot as plt
import logging
import numpy as np
import os
import sys

# Biopython modules
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.Data import SCOPData

# Non-standard modules
import contact_map as cm
import mutual_information as mut
import msa_caller as msa
import blast as blst
import plots

# Entrez inputs
Entrez.email = "ferran.muinos@gmail.com"
Entrez.tool = "cozmic.py"

# This program may well be called cozmic.py

if __name__ == '__main__':
	msg='Runs the main workflow'
	argparser = argparse.ArgumentParser(description=msg,
	            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = argparser.add_subparsers(help='Choose between modelize a pdb or from an existing file')
	# compulsory input: workflow type and target sequence/structure
# 	argparser.add_argument("workflow", help="""Type of workflow that will
# 	                       be conducted by the program: "real" will derive
# 	                       distances and contacts from the real PDB structure;
# 	                       "model" will do so with comparative modelling.""",
# 	                       choices=["real","model"])
	real = subparsers.add_parser('real', help='Make a structure of a real pdb')
	model = subparsers.add_parser('model', help='Modelize a pdb')
	
	real.add_argument("i", help="Target sequence id or filename")
	# argument options for cm module functions
	real.add_argument("-a", help="""Aa atom to calculate distance with:
                           CA (alpha carbon), CB (beta carbon), min (minimum 
                           distance between non-hydrogen atom pairs from each
                           residue).""",
                           choices=["CA", "CB", "min"],
                           default="min")
    # argument options for blst module functions
	real.add_argument("-blast", help="""Type of BLAST search to be 
	                       performed. Notice that blastp, blast and tblastn
	                       require a protein sequence, whereas blastx and 
	                       tblastx require a nucleotide sequence.""",
                           choices=["blastp","tblastn","blastx","tblastx"],
                           default="blastp")
                           # blastn does not return protein alignments
	choices_db = ["pdb", "swissprot", "refseq_protein", "pat",
                  "nr", "env_nr", "tsa_nr"]
	real.add_argument("-db", help="""Set the database of choice where
	                       BLAST will search.""",
                           choices=choices_db,
                           default="nr")
	real.add_argument("-s", help="""Set the maximum number hits 
	                       resulting from the BLAST search of homologs.""",
                           type=int,
                           default=200)
	real.add_argument("-f", help="""If present, then filter the BLAST 
	                       output by genus for attaining non-redundancy; 
	                       otherwise filter by species.""",
                           action='store_false', 
                           default=True)
    # argument options for msa and mut module functions
	real.add_argument("-g", help="""If present, then prune those 
                           columns of the MSA which have at least one gap.""",
                           action='store_true', 
                           default=False)
	# argument options for mut module functions
	real.add_argument("-b",
						   help="""Base of the logarithms used throughout
						   entropy and mutual information computations.""",
						   type=int,
                           default=20)
	real.add_argument("-low",
                           help="""Mininum entropy threshold allowed
                           for each column in the MSA. Columns below this
                           threshold are pruned for zMIc calculation.""",
                           type=float,
                           default=0.3)
	real.add_argument("-high",
                           help="""Maximum entropy threshold allowed 
                           for each column in the MSA. Columns above this
                           threshold are pruned for zMIc calculation.""",
                           type=float,
                           default=0.9)
	real.add_argument("-m",
                           help="""Method of choice for carrying out MSA.""",
                           default="clustalw",
                           choices=["clustalw", "muscle", "t_coffee"])
	args = argparser.parse_args()
	print(args)
	if args.i:
		# Retrieve the PDB structure, filter and get sequence
		parser = PDBParser(PERMISSIVE=1)
		if os.path.isfile(args.i):
			pdbpath = args.i
		else:
			pdbl = PDBList()
			try:
				pdbpath = pdbl.retrieve_pdb_file(args.i)
			except:
				sys.stderr.write("""exception: make sure your query format is correct.\n""")
		structure = parser.get_structure("cozmic_pdb_query", pdbpath)
		residues = cm.filter_residues(structure)
		s = ""
		for i in range(len(residues)): 
			s += SCOPData.protein_letters_3to1.get(residues[i].get_resname(), 'X') 
		seq = Seq(s,generic_protein)
		sys.stderr.write("Protein sequence:%s\n" %seq)
		# Compute distances and contact between residues
		dist_matrix = cm.calc_dist_matrix(residues, args.a)
		cont_matrix = cm.contact_map(dist_matrix, args.a)
		# Set up the files for calling BLAST
		blast_query_name = "cozmic_blast_query.fa"
		blast_out_name = "cozmic_blast.out"
		fd = open(blast_query_name,"w")
		fd.write(">%s\n%s" %("cozmic_blast_query",seq))
		fd.close()
		# Enter the query as first element in the output file
		file_out = open(blast_out_name, "w")
		record = SeqIO.read("cozmic_blast_query.fa", format="fasta")
		SeqIO.write(record, file_out, "fasta")
		# Call run_BLAST and write the output in the output file
		# So far it is not working
		blast_result = blst.run_BLAST(blast_query_name, args.blast, args.db, args.s)
		ides = blst.analyze_blast_result(blast_result, args.f)
		ids = list(blst.filter_ids(ides, "gi"))
		SeqIO.write(blst.retrive_sequence(ids), file_out, "fasta")
		file_out.close()
		# MSA: align the query to its homologs with the method of choice
		msa.call_msa_method(args.m, blast_out_name, "aligned.fa", "fasta")
		alignment = mut.AlignIO.read("aligned.fa", "fasta")
		# Prepare the alignment for MIc computations: 
		# prune high and low entropy columns
		edited = mut.prune_id_gaps(alignment, "cozmic_blast_query")
		gapped_list = []
		if args.g:
			gapped_list = mut.get_all_gaps(edited)
		edited = mut.prune(edited, gapped_list)
		(minlist, maxlist) = mut.get_extreme_columns(edited, args.low, args.high, args.b)
		edited = mut.prune(edited, minlist+maxlist)
		# compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
		MI_matrix = mut.mutual_info_matrix(edited, args.b)
		ncps_array = mut.NCPS_matrix(edited, args.b)
		MIc_matrix = MI_matrix - ncps_array
		zMIc_matrix = mut.standardise_matrix(MIc_matrix)
		# plot distance, contact, MIc Z-scores and its associated level matrix
		plots.plot_matrix_heatmap(dist_matrix,"Distance {}".format(args.a))
		plots.plot_matrix_binary(cont_matrix,"Contact {}".format(args.a))
		plots.plot_matrix_heatmap(zMIc_matrix,"zMIc")
		tmatrix = mut.get_level_matrix(zMIc_matrix,2)
		plots.plot_matrix_binary(tmatrix,"zMIc with L>2")
		# plot level-precision analysis and CM-distance analysis
		plots.precision_analysis(zMIc_matrix,cont_matrix,gapped_list,minlist,maxlist,0.0,3.0,60)
		# Leo's function here!
	elif args.pir:
		pass
		
	
	
	
