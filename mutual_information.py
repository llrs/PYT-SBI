#!/usr/bin/python3

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy
import copy
import math
import os

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
	matrix = numpy.empty([n,n], dtype=float, order='F')
	for j in range(n):
		for i in range(n):
			if j <= i: 
				matrix[i,j] = mutual_info(aln, i, j, base)
			else:
				matrix[i,j] = copy.copy(matrix[j,i])
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
	cps_matrix = numpy.empty([n,n], dtype=float, order='F')
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

	"""It returns the MIc matrix for aln"""
	
	return mutual_info_matrix(aln) - NCPS_matrix(aln)
			

if __name__ == "__main__":
	alignment = AlignIO.read("./data/alignlist.aln", "clustal")
	
#	alignment = MultipleSeqAlignment([
#            SeqRecord(Seq("MATGCTATTT", generic_protein), id="Alpha"),
#            SeqRecord(Seq("MAGTCTATGT", generic_protein), id="Beta"),
#            SeqRecord(Seq("MLTTCATTTG", generic_protein), id="Gamma"),
#            SeqRecord(Seq("MLGTCATTGG", generic_protein), id="Delta"),
#         ])
		
#	alignarray = numpy.array([list(rec) for rec in alignment], numpy.character)
#	print("Array shape %d by %d" %alignarray.shape)

#	alignarray = numpy.array([list(rec) for rec in edited], numpy.character)
#	print("Array shape %d by %d" %alignarray.shape)	

	edited = prune_first_gaps(alignment)
	print(edited)
	myarray = NCPS_matrix(edited,2)
	mi_matrix = mutual_info_matrix(edited,2)
	MIc_matrix = mi_matrix-myarray
	print(mi_matrix)
	print(myarray)
	print(MIc_matrix)
	
	# plot 1
	plt.imshow(mi_matrix, interpolation='none')
	heatmap = plt.pcolormesh(mi_matrix)
	plt.title('Fancy MI Heatmap')
	legend = plt.colorbar(heatmap)
	legend.set_label("MI")
	plt.savefig('fancy_heatmap_file1.png')
	fig = plt.figure()
	fig.show()

	# plot 2
	plt.imshow(myarray, interpolation='none')
	heatmap = plt.pcolormesh(myarray)
	plt.title('Fancy NCPS Heatmap')
	legend = plt.colorbar(heatmap)
	legend.set_label("NCPS")
	plt.savefig('fancy_heatmap_file2.png')
	fig = plt.figure()
	fig.show()
	
	# plot 3
	plt.imshow(MIc_matrix, interpolation='none')
	heatmap = plt.pcolormesh(MIc_matrix)
	plt.title('Fancy MIc Heatmap')
	legend = plt.colorbar(heatmap)
	legend.set_label("MIc")
	plt.savefig('fancy_heatmap_file3.png')
	fig = plt.figure()
	fig.show()
