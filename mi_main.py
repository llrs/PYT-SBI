#!/usr/bin/python3
'''
Created on Mar 7, 2016
@author: Lluís, Leo, Ferran
'''
import argparse
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
	
#   base = mut.os.path.basename(args.f)
#   name_f = mut.os.path.splitext(base)[0]
    alignment = mut.AlignIO.read(args.f, "clustal")
    
    edited = mut.prune_first_gaps(alignment)
    gapped_list = mut.get_all_gaps(edited)
    if args.g != '0':
        edited = mut.prune(edited,gapped_list)
    (minlist, maxlist) = mut.get_extreme_columns(edited,args.low,args.high,args.b)
    edited = mut.prune(edited,minlist+maxlist)
	
    MI_matrix = mut.mutual_info_matrix(edited,args.b)
    ncps_array = mut.NCPS_matrix(edited,args.b)
    MIc_matrix = MI_matrix - ncps_array
    zMIc_matrix = mut.standardise_matrix(MIc_matrix)
    tmatrix = mut.get_level_matrix(zMIc_matrix,args.l)

    mut.plot_matrix_heatmap(zMIc_matrix,"zMIc")
    mut.plot_matrix_binary(tmatrix,"zMIc_level")
	
    cm_residue_positions = mut.retrieve_residue_positions(tmatrix,gapped_list,minlist+maxlist)
    print(cm_residue_positions)
