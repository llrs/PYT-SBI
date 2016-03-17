#!/usr/bin/python3
"""
Deals with computations related to mutual information, phylogenetic
signals and structural/functional correlation signals from data 
encoded in multiple sequence alignments.
Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import argparse
import matplotlib.pyplot as plt
import numpy as np
import copy
import math
import os

# Biopython modules
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Non-standard modules
import plots

def prune(aln, listcol):
    """It prunes the colunms enclosed in listcol from the
    input MultipleSeqAlignment object aln."""
    listcol = sorted(listcol)
    for i in range(len(listcol)):
        aln = aln[:, : listcol[i]-i] + aln[:, listcol[i]-i+1:]
    return aln


def prune_id_gaps(aln, identifier):
    """Prunes the gaps in aln record labeled with identifier"""
    for i in range(len(aln)):
        if aln[i].id == identifier:
            row = i
    listcol = []
    for j in range(aln.get_alignment_length()):
        if aln[row].seq[j] == '-':
            listcol.append(j)
    return prune(aln, listcol)


def prune_first_gaps(aln):
    """It returns a MultipleSeqAlignment object where the columns of
    the input with a gap in the first row have been removed."""
    listcol = []
    for i in range(aln.get_alignment_length()):
        if aln[0].seq[i] == '-':
            listcol.append(i)
    return prune(aln, listcol)


def get_all_gaps(aln):
    """It returns a list of columns of aln with gaps at any position"""
    listcol = []
    for j in range(aln.get_alignment_length()):
        for i in range(len(aln)):
            if aln[i].seq[j] == '-':
                listcol.append(j)
                break
    return listcol


def get_level_matrix(matrix, level):
    """It returns a matrix with binary values representing the cells
    going higher the threshold level in input matrix"""
    (n1, n2) = matrix.shape
    out_matrix = np.empty([n1, n2], dtype=float, order='F')
    for i in range(n1):
        for j in range(n2):
            if i == j:
                out_matrix[i, j] = 0
            elif matrix[i, j] >= level:
                out_matrix[i, j] = 1
            else:
                out_matrix[i, j] = 0
    return out_matrix


def column_frequencies(aln, col):
    """ It returns a dictionary with the frequency for each residue
    in column col of MultipleSeqAlignment aln"""
    alphabet = set('ACDEFGHIKLMNPQRSTVWY-')
    freq = dict.fromkeys(alphabet, 0)
    column = aln[:, col]
    for char in column:
        try:
            freq[char] += 1
        except KeyError:
            pass
    for key in freq:
        freq[key] = freq[key]/len(column)
    return freq


def joint_column_frequencies(aln, col1, col2):
    """It returns a dictionary with the joint frequency for each
    ordered pair (tuple) of residues in columns indexed col1 and col2,
    respectively, of the MultipleSeqAlignment aln"""
    alphabet = set('ACDEFGHIKLMNPQRSTVWY-')
    index = set()
    for i in alphabet:
        for j in alphabet:
            index.add((i, j))
    freq = dict.fromkeys(index, 0)
    columnx = aln[:, col1]
    columny = aln[:, col2]
    for i in range(len(aln)):
        try:
            freq[(columnx[i], columny[i])] += 1
        except KeyError:
            pass
    for key in freq:
        freq[key] = freq[key]/len(aln)
    return freq


def entropy(aln, col, base):
    """It returns the entropy for column col in the MultipleSeqAlignment
    object aln"""
    freq = column_frequencies(aln, col)
    entropy = 0
    for key in freq:
        if freq[key] != 0:
            entropy -= freq[key]*math.log(freq[key], base)
    return entropy


def get_extreme_columns(aln, entmin, entmax, base):
    """It returns a list of columns of aln with entropy below the threshold"""
    minlist = []
    maxlist = []
    for i in range(aln.get_alignment_length()):
        if entropy(aln, i, base) <= entmin:
            minlist.append(i)
        if entropy(aln, i, base) >= entmax:
            maxlist.append(i)
    return (minlist, maxlist)


def joint_entropy(aln, col1, col2, base):
    """It returns the joint entropy for col1 and col2 in the
    MultipleSeqAlignment object aln"""
    joint_freq = joint_column_frequencies(aln, col1, col2)
    jentropy = 0
    for key in joint_freq:
        if joint_freq[key] != 0:
            jentropy -= joint_freq[key]*math.log(joint_freq[key], base)
    return jentropy


def mutual_info(aln, col1, col2, base=20):
    """ It returns the MI of a pair of columns"""
    entropy1 = entropy(aln, col1, base)
    entropy2 = entropy(aln, col2, base)
    return entropy1 + entropy2 - joint_entropy(aln, col1, col2, base)


def mutual_info_matrix(aln, base=20):
    """It returns a square matrix of size len(aln) with MI
    values for each pair of positions in the MultipleSeqAlignment
    object aln"""
    n = aln.get_alignment_length()
    matrix = np.empty([n, n], dtype=float, order='F')
    for j in range(n):
        for i in range(n):
            if j <= i:
                matrix[i, j] = mutual_info(aln, i, j, base)
            else:
                matrix[i, j] = copy.copy(matrix[j, i])
    return matrix

def matrix_hits(binary_matrix):
	"""Number of True cells above the diagonal in a binary matrix"""
	count = 0
	(n1,n2) = binary_matrix.shape
	for i in range(n1):
		for j in range(n2):
			if j > i:
				count += binary_matrix[i,j]
	return count

def standardise_matrix(mat):
    """It resturns a matrix with Z-score values of mat, using the mean and
    standard deviation estimators over all the entries of mat"""
    myarray = []
    (n1, n2) = mat.shape
    for i in range(n1):
        for j in range(n2):
            if i != j:
                myarray.append(mat[i, j])
    mymean = np.mean(myarray)
    mystd = np.std(myarray)
    matrix = np.empty([n1, n2], dtype=float, order='F')
    for i in range(n1):
        for j in range(n2):
            if i == j:
                matrix[i, j] = 0
            else:
                matrix[i, j] = (mat[i, j]-mymean)/mystd
    return matrix


def CPS(aln, col1, col2, base=20):
    """It returns the pairwise co-evolutionary pattern similarity
    for a pair of columns col1 and col2 of the MultipleSeqAlignment
    object aln"""
    n = aln.get_alignment_length()
    m = mutual_info_matrix(aln, base)
    cps = 0
    for k in range(n):
        if k != col1 and k != col2:
            cps += m[col1, k]*m[k, col2]
    cps = cps/(n-2)
    return cps


def NCPS_matrix(aln, base=20):
    """It returns the matrix of pairwise normalised co-evolutionary
    pattern similarities of the MuplipleSeqAlignment aln"""
    n = aln.get_alignment_length()
    m = mutual_info_matrix(aln, base)
    cps_matrix = np.empty([n, n], dtype=float, order='F')
    den = 0
    for j in range(n):
        for i in range(n):
            cps = 0
            if j <= i:
                for k in range(n):
                    if k != i and k != j:
                        cps += m[i, k]*m[k, j]
                cps_matrix[i, j] = cps / (n-2)
            else:
                cps_matrix[i, j] = copy.copy(cps_matrix[j, i])
            den += cps_matrix[i, j]
    den = den/(n*(n-1))
    den = math.sqrt(den)
    return cps_matrix * 1/den


def mutual_info_c(aln):
    """It returns the MIc and its standardised version from aln"""
    mic = mutual_info_matrix(aln) - NCPS_matrix(aln)
    return (mic, standardise_matrix(mic))


def plot_matrix_heatmap(matrix, keyword):
    """Plots a matrix of values as a heatmap"""

    imgplot = plt.imshow(matrix, cmap='Blues', interpolation='none')
    # heatmap = plt.pcolormesh(matrix, cmap='Blues')
    plt.title('{} heatmap'.format(keyword))
    legend = plt.colorbar(imgplot)
    legend.set_label("{}".format(keyword))
    plt.savefig('{}_heatmap.png'.format(keyword), format="png")
    fig = plt.figure()
    fig.show()


def plot_matrix_binary(matrix, keyword):
    """Plots a matrix with binary values"""

    imgplot = plt.imshow(matrix, cmap = 'Greys', interpolation='none')
    plt.title('{}'.format(keyword))
    plt.savefig('{}_heatmap.png'.format(keyword), format="png")
    fig = plt.figure()
    fig.show()


def reconstruct_position(pos, deleted_pos):
    """Returns the original position of a column, currently in position
    pos, in an MSA that has undergone prunning in columns at the
    positions comprised in the list deleted_pos."""
    diff = 0
    deleted_pos = sorted(deleted_pos)
    for i in range(len(deleted_pos)):
        if deleted_pos[i] <= pos + diff:
            diff += 1
        else:
            break
    return pos + diff


def retrieve_residue_positions(binary_matrix, gap_list, extreme_list):
    """Retrieves the original coordinates of the set of pairs of residues
    corresponding to True spots in the binary matrix; only those above
    the diagonal suffice"""
    set_pairs = set()
    (n1, n2) = binary_matrix.shape
    for i in range(n1):
        i_reconstruct = reconstruct_position(i, extreme_list)
        for j in range(n2):
            if j > i and binary_matrix[i, j] == 1:
                j_reconstruct = reconstruct_position(j, extreme_list)
                res1 = reconstruct_position(i_reconstruct, gap_list)
                res2 = reconstruct_position(j_reconstruct, gap_list)
                set_pairs.add(tuple(sorted([res1, res2])))
    return set_pairs


def retrieve_all_positions(matrix, gap_list, extreme_list):
    """Retrieves a dictionary with keys given by the original coordinates
    of the set of pairs of residues corresponding to all spots in matrix
    above the diagonal, and values corresponding to the values in matrix."""
    dict_pairs = {}
    (n1, n2) = matrix.shape
    for i in range(n1):
        i_reconstruct = reconstruct_position(i, extreme_list)
        for j in range(n2):
            if j > i:
                j_reconstruct = reconstruct_position(j, extreme_list)
                res1 = reconstruct_position(i_reconstruct, gap_list)
                res2 = reconstruct_position(j_reconstruct, gap_list)
                value = matrix[i, j]
                dict_pairs[tuple(sorted([res1, res2]))] = value
    return dict_pairs

if __name__ == "__main__":
	msg='Runs the computations related to zMIc'
	argparser = argparse.ArgumentParser(description=msg,
	            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# compulsory input
	argparser.add_argument("i", help="MSA file")
	argparser.add_argument("id", help="Target identifier")
	# optional input
	argparser.add_argument("-o",
						   help="""Output file to keep the CM residue
						   pair candidates.""",
						   default="CM_residue_pairs.out")
	argparser.add_argument("-L",
						   help="""Threshold level to select candidates.""",
						   type=float,
						   default=2.0)
	argparser.add_argument("-b",
						   help="""Base of the logarithms used throughout
						   entropy and mutual information computations.""",
						   type=int,
                           default=20)
	argparser.add_argument("-g", help="""If present, then prune those 
                           columns of the MSA which have at least one gap.""",
                           action='store_true', 
                           default=False)
	argparser.add_argument("-low",
                           help="""Threshold of mininum entropy allowed
                           for each column in a MSA. Columns below this
                           threshold are pruned.""",
                           type=float,
                           default=0.3)
	argparser.add_argument("-high",
                           help="""Threshold of maximum entropy allowed 
                           for each column in a MSA. Columns above this
                           threshold are pruned.""",
                           type=float,
                           default=0.9)
	args = argparser.parse_args()
	# Read MSA
	alignment = AlignIO.read(args.i, "clustal")
	# Prepare the alignment for MIc computations: prune high and low entropy columns
	edited = prune_id_gaps(alignment, args.id)
	gapped_list = []
	if args.g:
		gapped_list = get_all_gaps(edited)
	edited = prune(edited, gapped_list)
	(minlist, maxlist) = get_extreme_columns(edited, args.low, args.high, args.b)
	edited = prune(edited, minlist+maxlist)
	# compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
	MI_matrix = mutual_info_matrix(edited, args.b)
	ncps_array = NCPS_matrix(edited, args.b)
	MIc_matrix = MI_matrix - ncps_array
	zMIc_matrix = standardise_matrix(MIc_matrix)
	# plot MIc Z-scores and its associated level matrix
	plots.plot_matrix_heatmap(zMIc_matrix,"zMIc")
	tmatrix = get_level_matrix(zMIc_matrix,args.L)
	plots.plot_matrix_binary(tmatrix,"zMIc>{}".format(args.L))
	# Retrieve CM residue pairs in their original coordinates
	cm_residue_pairs = retrieve_residue_positions(tmatrix,gapped_list,minlist+maxlist)
	fd = open(args.o, "w")
	for residue_pair in sorted(cm_residue_pairs):
		fd.write("%d %d\n" %residue_pair)
	fd.close()
