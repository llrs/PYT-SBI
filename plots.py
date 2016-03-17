#!/usr/bin/python3
"""
Deals with plots and graphical outputs.
Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import matplotlib.pyplot as plt
import numpy as np
import logging

# non-standard modules
import contact_map as cm
import mutual_information as mut

def plot_distance(distances, name_file, option):
    """Plots the distances between the residues."""
    logging.info("Plotting the distance map for {}".format(name_f))
    plt.imshow(distances, interpolation='none')
    heatmap = plt.pcolormesh(distances)
    plt.title('Distances of the file {}'.format(name_file))
    legend = plt.colorbar(heatmap)
    legend.set_label("Angstroms")
    plt.savefig('distance_map_{}_{}.png'.format(name_file, option),
                format="png")

def plot_contacts(contacts, name_file, option):
    """Plots the contact map between residues"""
    logging.info("Plotting the contact map for {}".format(name_file))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle("Contact between residues of {}".format(name_file))
    ax.imshow(contacts, aspect='auto',
              cmap=plt.cm.gray, interpolation='nearest')
    plt.savefig("contact_map_{}_{}.png".format(name_file, option),
                format="png")

def plot_matrix_heatmap(matrix, keyword):
	"""Plots a matrix of values as a heatmap"""
	imgplot = plt.imshow(matrix, cmap='Blues', interpolation='none')
	plt.title('{} heatmap'.format(keyword))
	legend = plt.colorbar(imgplot)
	legend.set_label("{}".format(keyword))
	plt.savefig('{}_heatmap.png'.format(keyword), format="png")
	fig = plt.figure()
	fig.show()

def plot_matrix_binary(matrix, keyword):
	"""Plots a matrix with binary values: black and white"""
	imgplot=plt.imshow(matrix, cmap = 'Greys', interpolation='none')
	plt.title('{}'.format(keyword))
	plt.savefig('{}_heatmap.png'.format(keyword), format="png")
	fig = plt.figure()
	fig.show()

def plot_twin_curves(cutoff_list, hit_list, precision_list):
	"""
	Plots the cutoff_list x hit_list and cutoff_list x precision_list 
	curves in the same plot with two different yaxis labels.
	"""
	fig, ax1 = plt.subplots()
	ax1.plot(cutoff_list, hit_list,'ob-')
	ax1.set_xlabel('Cutoff level L')
	# Make the y-axis label match the line color.
	ax1.set_ylabel('Predicted CM residue pairs', color='b')
	# Make the second plot with another axis
	ax2 = ax1.twinx()
	ax2.plot(cutoff_list, precision_list,'or-')
	ax2.set_ylabel('Precision', color='r')
	plt.show()

def precision_analysis(zMIc_matrix,cont_matrix,gapped_list,minlist,maxlist,l=0.0,h=3.0,num=60):
	"""
	Given zMIc and contacts, computes and associates the precision
	with the number of CM predictions with the threshold level 
	"""
	cutoff_list = []
	hit_list = []
	precision_list = []
	for cutoff in np.linspace(l,h,num):
		cutoff_list.append(cutoff)
		tmatrix = mut.get_level_matrix(zMIc_matrix,cutoff)
		hits = mut.matrix_hits(tmatrix)
		hit_list.append(hits)
		cm_residue_pairs = mut.retrieve_residue_positions(tmatrix,gapped_list,minlist+maxlist)
		print(cm_residue_pairs)
		count = 0
		for rp in cm_residue_pairs:
			count += cont_matrix[rp[0],rp[1]]
		precision_list.append(count/hits)
	plot_twin_curves(cutoff_list, hit_list, precision_list)

if __name__ == "__main__":
	pass
	
