#!/usr/bin/python3
'''
Provides the workflow from an input MSA all the way to a standardised 
form of MIc matrix and CM residue pair candidates list thereof
Created on Mar 12, 2016
@author: Leo, Ferran, Lluís
'''
# standard modules
import argparse
import os

# non-standard own module
import mutual_information as mut

if __name__ == '__main__':
	msg = 'Extracts correlated mutation signals from MSA'
	argparser = argparse.ArgumentParser(description=msg)
	argparser.add_argument("-f", help="MSA to analyze.")
	argparser.add_argument("-g", 
                           help="Prunning of column gaps.",
                           default=0)
	argparser.add_argument("-b", 
						   help="Base of the logarithms.",
                           default=20)
	argparser.add_argument("-low", 
                           help="Min entropy threshold.",
                           default=0.3)
	argparser.add_argument("-high", 
                           help="Max entropy threshold.",
                           default=0.9)
	argparser.add_argument("-l", 
                           help="Z-score CM level in standard deviations.",
                           default=2)
    
	args = argparser.parse_args()
	
	base = os.path.basename(args.f)
	name_f = os.path.splitext(base)[0]
	alignment = mut.AlignIO.read(args.f, "clustal")
    
	edited = mut.prune_first_gaps(alignment)
	gapped_list = mut.get_all_gaps(edited)
	if args.g != '0':
		edited = mut.prune(edited,gapped_list)
	(minlist, maxlist) = mut.get_extreme_columns(edited,args.low,args.high,args.b)
	edited = mut.prune(edited,minlist+maxlist)
	
	# compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
	MI_matrix = mut.mutual_info_matrix(edited,args.b)
	ncps_array = mut.NCPS_matrix(edited,args.b)
	MIc_matrix = MI_matrix - ncps_array
	zMIc_matrix = mut.standardise_matrix(MIc_matrix)
	tmatrix = mut.get_level_matrix(zMIc_matrix,args.l)
	
	# plot MIc Z-scores and its associated level matrix
	mut.plot_matrix_heatmap(zMIc_matrix,"zMIc")
	mut.plot_matrix_binary(tmatrix,"zMIc_level")
	
	# return CM residue pair candidates in lexicographic order
	cm_residue_pairs = mut.retrieve_residue_positions(tmatrix,gapped_list,minlist+maxlist)
	fd = open(name_f + ".out", "w")
	for residue_pair in sorted(cm_residue_pairs):
		fd.write("%d %d\n" %residue_pair)
	fd.close()
