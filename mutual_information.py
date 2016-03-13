#!/usr/bin/python3
"""
Deals with computations related to mutual information, phylogenetic
signals and structural/functional correlation signals from data 
encoded in multiple sequence alignments.

Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import matplotlib.pyplot as plt
import numpy as np
import copy
import math
import os

# Biopython modules
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def prune(aln, listcol):
	""" 
	It prunes the colunms enclosed in listcol from the 
	input MultipleSeqAlignment object aln.
	"""
	listcol = sorted(listcol)
	for i in range(len(listcol)):
		aln = aln[:,:listcol[i]-i] + aln[:,listcol[i]-i+1:]
	return aln

def prune_first_gaps(aln):
	""" 
	It returns a MultipleSeqAlignment object where the columns of 
	the input with a gap in the first row have been removed.
	"""
	listcol = []
	for i in range(aln.get_alignment_length()):
		if aln[0].seq[i] == '-':
			listcol.append(i)
	return prune(aln,listcol)
	
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
	""" 
	It returns a matrix with binary values representing the cells
	going higher the threshold level in input matrix
	""" 
	(n1, n2) = matrix.shape
	out_matrix = np.empty([n1,n2], dtype=float, order='F')
	for i in range(n1):
		for j in range(n2):
			if i == j:
				out_matrix[i,j] = 0
			elif matrix[i,j] >= level:
				out_matrix[i,j] = 1
			else:
				out_matrix[i,j] = 0
	return out_matrix

def column_frequencies(aln, col):
	""" 
	It returns a dictionary with the frequency for each residue 
	in column col of MultipleSeqAlignment aln
	"""
	alphabet = set('ACDEFGHIKLMNPQRSTVWY-')
	freq = dict.fromkeys(alphabet,0)
	column = aln[:,col]
	for char in column:
		freq[char] += 1
	for key in freq:
		freq[key] = freq[key]/len(column)
	return freq

def joint_column_frequencies(aln, col1, col2):
	""" 
	It returns a dictionary with the joint frequency for each 
	ordered pair (tuple) of residues in columns indexed col1 and col2, 
	respectively, of the MultipleSeqAlignment aln
	"""
	alphabet = set('ACDEFGHIKLMNPQRSTVWY-')
	index = set()
	for i in alphabet:
		for j in alphabet:
			index.add((i,j))
	freq = dict.fromkeys(index,0)
	columnx = aln[:,col1]
	columny = aln[:,col2]
	for i in range(len(aln)):
		freq[(columnx[i],columny[i])] += 1
	for key in freq:
		freq[key] = freq[key]/len(aln)
	return freq

def entropy(aln, col, base):
	""" 
	It returns the entropy for column col in the MultipleSeqAlignment
	object aln 
	"""
	freq = column_frequencies(aln, col)
	entropy = 0
	for key in freq:
		if freq[key] != 0:
			entropy -= freq[key]*math.log(freq[key],base)
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
	""" 
	It returns the joint entropy for col1 and col2 in the 
	MultipleSeqAlignment object aln """
	joint_freq = joint_column_frequencies(aln, col1, col2)
	jentropy = 0
	for key in joint_freq:
		if joint_freq[key] != 0:
			jentropy -= joint_freq[key]*math.log(joint_freq[key],base)
	return jentropy
	
def mutual_info(aln, col1, col2, base=20):
	""" It returns the MI of a pair of columns"""
	return entropy(aln, col1, base) + entropy(aln, col2, base) - joint_entropy(aln, col1, col2, base)
	
def mutual_info_matrix(aln, base=20):
	""" 
	It returns a square matrix of size len(aln) with MI 
	values for each pair of positions in the MultipleSeqAlignment
	object aln
	"""
	n = aln.get_alignment_length()
	matrix = np.empty([n,n], dtype=float, order='F')
	for j in range(n):
		for i in range(n):
			if j <= i: 
				matrix[i,j] = mutual_info(aln, i, j, base)
			else:
				matrix[i,j] = copy.copy(matrix[j,i])
	return matrix

def standardise_matrix(mat):
	"""
	It resturns a matrix with Z-score values of mat, using the mean and 
	standard deviation estimators over all the entries of mat 
	"""
	myarray = []
	(n1,n2) = mat.shape
	for i in range(n1):
		for j in range(n2):
			if i != j:
				myarray.append(mat[i,j])
	mymean = np.mean(myarray)
	mystd = np.std(myarray)
	matrix = np.empty([n1,n2], dtype=float, order='F')
	for i in range(n1):
		for j in range(n2):
			if i == j:
				matrix[i,j] = 0
			else:
				matrix[i,j] = (mat[i,j]-mymean)/mystd
	return matrix

def CPS(aln, col1, col2, base=20):
	""" 
	It returns the pairwise co-evolutionary pattern similarity
	for a pair of columns col1 and col2 of the MultipleSeqAlignment
	object aln
	"""
	n = aln.get_alignment_length()
	m = mutual_info_matrix(aln,base)
	cps = 0
	for k in range(n):
		if k != col1 and k != col2:
			cps += m[col1,k]*m[k,col2]
	cps = cps/(n-2)
	return cps
	
def NCPS_matrix(aln, base=20):
	""" 
	It returns the matrix of pairwise normalised co-evolutionary 
	pattern similarities of the MuplipleSeqAlignment aln
	"""
	n = aln.get_alignment_length()
	m = mutual_info_matrix(aln,base)
	cps_matrix = np.empty([n,n], dtype=float, order='F')
	den = 0
	for j in range(n):
		for i in range(n):
			cps = 0
			if j <= i:
				for k in range(n):
					if k != i and k != j: 
						cps += m[i,k]*m[k,j]
				cps_matrix[i,j] = cps / (n-2)
			else:
				cps_matrix[i,j] = copy.copy(cps_matrix[j,i])
			den += cps_matrix[i,j]
	den = den/(n*(n-1))
	den = math.sqrt(den)
	return cps_matrix * 1/den

def mutual_info_c(aln):
	"""It returns the MIc and its standardised version from aln"""
	mic = mutual_info_matrix(aln) - NCPS_matrix(aln)
	return (mic, standardise_matrix(mic))


def plot_matrix_heatmap(matrix, keyword):
	"""Plots a matrix of values as a heatmap"""
	
	plt.imshow(matrix, interpolation='none')
	heatmap = plt.pcolor(matrix,cmap=plt.cm.binary)
	plt.title('{} heatmap'.format(keyword))
	legend = plt.colorbar(heatmap)
	legend.set_label("{}".format(keyword))
	plt.savefig('{}_heatmap.png'.format(keyword), format="png")
	fig = plt.figure()
	fig.show()

def plot_matrix_binary(matrix, keyword):
	"""Plots a matrix with binary values"""
	
	plt.imshow(matrix, interpolation='none')
	heatmap = plt.pcolor(matrix,cmap=plt.cm.binary)
	plt.title('{} heatmap'.format(keyword))
	plt.savefig('{}_heatmap.png'.format(keyword), format="png")
	fig = plt.figure()
	fig.show()

def reconstruct_position(pos, deleted_pos):
	"""
	Returns the original position of a column, currently in position 
	pos, in an MSA that has undergone prunning in columns at the 
	positions comprised in the list	deleted_pos.
	"""
	diff = 0
	deleted_pos = sorted(deleted_pos)
	for i in range(len(deleted_pos)):
		if deleted_pos[i] <= pos + diff:
			diff += 1
		else:
			break
	return pos + diff
			
def retrieve_residue_positions(binary_matrix, gap_list, extreme_list):
	"""
	Retrieves the set of pairs of residues showing signals of 
	correlated mutations from the binary level matrix
	"""
	set_pairs = set()
	(n1,n2) = binary_matrix.shape
	for i in range(n1):
		for j in range(n2):
			if binary_matrix[i,j] == 1:
				res1 = reconstruct_position(reconstruct_position(i,extreme_list),gap_list)
				res2 = reconstruct_position(reconstruct_position(j,extreme_list),gap_list)
				set_pairs.add(frozenset([res1,res2]))
	return set_pairs
				
if __name__ == "__main__":
#	alignment = AlignIO.read("./data/alignlist.aln", "clustal")
	
	alignment = MultipleSeqAlignment([
            SeqRecord(Seq("MATTGCTATTT", generic_protein), id="Alpha"),
            SeqRecord(Seq("MAGT-CTATGT", generic_protein), id="Beta"),
            SeqRecord(Seq("MLTTCA-TTTG", generic_protein), id="Gamma"),
            SeqRecord(Seq("-LGTCATTG-G", generic_protein), id="Delta"),
         ])

	edited = prune(alignment, [0,4,6,9])
	print(edited)
	print(reconstruct_position(6,[]))	
	
#	alignarray = np.array([list(rec) for rec in edited], np.character)
#	print("Array shape %d by %d" %alignarray.shape)	

