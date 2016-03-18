#!/usr/bin/python3
# encoding: utf-8
"""
Deals with computations related to mutual information, phylogenetic
signals and structural/functional correlation signals from data
encoded in multiple sequence alignments.
Created on Mar 7, 2016
@author: Leo, Lluís, Ferran
"""
# standard modules
import argparse
import copy
import math
import logging

import numpy as np

# Biopython modules
from Bio import AlignIO

# Non-standard modules
import plots


def prune(aln, listcol):
    """Removes columns from the multiple sequence alignment (MSA).

    aln: MultipleSeqAlignment object
    list: list of indexes of columns of aln
    It removes the columns listed in listcol from aln."""

    logging.debug("Pruning alignment.")

    listcol = sorted(listcol)
    for i in range(len(listcol)):
        aln = aln[:, : listcol[i]-i] + aln[:, listcol[i]-i+1:]
    return aln


def prune_id_gaps(aln, identifier):
    """Removes columns which have gaps in a specific record.

    aln: MultipleSeqAlignment object
    identifer: string of the record's identifier
    It removes the columns of aln which have gaps in the record
    labeled with identifier."""

    logging.debug("Pruning gaps by id.")

    for i in range(len(aln)):
        if aln[i].id == identifier:
            row = i
    listcol = []
    for j in range(aln.get_alignment_length()):
        if aln[row].seq[j] == '-':
            listcol.append(j)
    return prune(aln, listcol)


def get_all_gaps(aln):
    """Gets columns with at least one gap."""

    logging.info("Identifying positions with gaps.")

    listcol = []
    for j in range(aln.get_alignment_length()):
        for i in range(len(aln)):
            if aln[i].seq[j] == '-':
                listcol.append(j)
                break
    return listcol


def get_level_matrix(matrix, level):
    """Returns a binary matrix with positions exceeding a threshold.

    matrix = numpy array object
    level = floating number
    The matrix it returns has 1 in the positions where matrix
    has values above level and 0 elsewhere."""

    logging.info("Selecting the amino acids contacts.")

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
    """ Computes the residue frequencies of an MSA column.

    aln: MultipleSeqAlignment object
    col: integer index of a column
    It returns a dictionary with the frequency (value) for each residue
    in column col of the MultipleSeqAlignment aln."""

    logging.info("Calculating the frecuency for the alignment.")

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
    """Computes the joint frequencies of pairs of MSA columns.

    aln: MultipleSeqAlignment object
    col1, col2: integer indexes of columns
    It returns a dictionary with the joint frequency (value) for each
    sorted tuple (key) of residues in columns indexed col1 and col2
    of the MultipleSeqAlignment aln."""

    msg_log = "Calculates the joint frequency of columns {} and {}"
    logging.info(msg_log.format(col1, col2))

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
    """Computes the entropy of an MSA column.

    aln: MultipleSeqAlignment object
    col: integer index of a column
    It returns the entropy for column col in aln."""

    logging.info("Calculating the entropy at column {}".format(col))

    freq = column_frequencies(aln, col)
    entropy = 0
    for key in freq:
        if freq[key] != 0:
            entropy -= freq[key]*math.log(freq[key], base)
    return entropy


def joint_entropy(aln, col1, col2, base):
    """Computes the joint entropy of a pair of MSA columns.

    aln: MultipleSeqAlignment object
    col1, col2: integer indexes of columns
    base:  base of the logarithm
    It returns the joint entropy for col1 and col2 in aln."""
    logging.info("Calculates joint entropy of {} and {}".format(col1, col2))
    joint_freq = joint_column_frequencies(aln, col1, col2)
    jentropy = 0
    for key in joint_freq:
        if joint_freq[key] != 0:
            jentropy -= joint_freq[key]*math.log(joint_freq[key], base)
    return jentropy


def get_extreme_columns(aln, entmin, entmax, base):
    """Gets columns with extreme entropy values.

    aln: MultipleSeqAlignment object
    entmin, entmax: min (resp. max) entropy thresholds
    base: base of the logarithm
    It returns a list of column indexes with entropy below (resp. above)
    the thresholds"""

    msg_log = "Calculating the columns with entropy below the threshold."
    logging.info(msg_log)

    minlist = []
    maxlist = []
    for i in range(aln.get_alignment_length()):
        if entropy(aln, i, base) <= entmin:
            minlist.append(i)
        if entropy(aln, i, base) >= entmax:
            maxlist.append(i)
    return (minlist, maxlist)


def mutual_info(aln, col1, col2, base=20):
    """Gives the mutual information (MI) of a pair of MSA columns.

    aln: MultipleSeqAlignment object
    col1, col2: integer indexes of columns
    base:  base of the logarithm"""

    msg_log = "Calculates the mutual information of {} and {} columns"
    logging.info(msg_log.format(col1, col2))

    entropy1 = entropy(aln, col1, base)
    entropy2 = entropy(aln, col2, base)
    return entropy1 + entropy2 - joint_entropy(aln, col1, col2, base)


def mutual_info_matrix(aln, base=20):
    """Returns the mutual information (MI) matrix of an MSA.

    aln: MultipleSeqAlignment object
    base: base of the logarithm
    It returns a square matrix of size len(aln) with MI
    values for each pair of positions in the MultipleSeqAlignment
    object aln"""

    msg_log = "Calculates the matrix of the mutual information"
    logging.info(msg_log)

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
    """Gets the number of cells with value 1 in matrix."""

    logging.info("Counting how many contacts are predicted.")

    count = 0
    (n1, n2) = binary_matrix.shape
    for i in range(n1):
        for j in range(n2):
            if j > i:
                count += binary_matrix[i, j]
    return count


def standardise_matrix(mat):
    """Translates the values of matrix into Z-scores.

    mat: numpy array of floats
    It returns a matrix with Z-score values of mat, using the mean and
    standard deviation estimators over all the entries of mat, excluding
    the values in the diagonal."""

    logging.info("Standardising the MI with Z-score.")

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
    """Computes the CPS of a pair of columns in an MSA.

    aln: MultipleSeqAlignment object
    col1, col2: integer indexes of columns
    base:  base of the logarithm
    It returns the pairwise co-evolutionary pattern similarity
    for a pair of columns col1 and col2 of aln"""

    logging.info("Calculating the co-evolutionary pattern similarity.")

    n = aln.get_alignment_length()

    if n <= 1:
        raise ValueError("Too few sequence in the alignment provided.")

    m = mutual_info_matrix(aln, base)
    cps = 0
    for k in range(n):
        if k != col1 and k != col2:
            cps += m[col1, k]*m[k, col2]
    cps = cps/(n-2)
    return cps


def NCPS_matrix(aln, base=20):
    """Computes the NCPS matrix of a MSA.

    aln: MultipleSeqAlignment object
    base:  base of the logarithm
    It returns the matrix of pairwise normalised co-evolutionary
    pattern similarities of pairs of columns in aln"""

    logging.info("Normalizing co-evolutionary pattern similarities")

    n = aln.get_alignment_length()

    if n <= 1:
        raise ValueError("Too few sequence in the alignment provided. Length of the aligment {}".format(n))

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
    """It returns the MIc and its standardised version from aln."""

    logging.info("Calculating the MIc and standard MI.")

    mic = mutual_info_matrix(aln) - NCPS_matrix(aln)
    return (mic, standardise_matrix(mic))


def reconstruct_position(pos, deleted_pos):
    """Computes the prior position of one column in an MSA.

    pos: index of current position in an MSA
    deleted_positions: list with column indexes that were removed
    It returns the position of a column in an MSA that has undergone
    deletion of several columns; the list deleted_pos contains the
    indexes of the columns that have been deleted, prior to deletion."""

    logging.info("Reconstructing original positions of the matrix.")

    diff = 0
    deleted_pos = sorted(deleted_pos)
    for i in range(len(deleted_pos)):
        if deleted_pos[i] <= pos + diff:
            diff += 1
        else:
            break
    return pos + diff


def retrieve_residue_positions(binary_matrix, gap_list, extreme_list):
    """Computes the prior indexes of column pairs.

    pos: index of current position in an MSA
    gap_list, extreme_list: lists of indexes of column that were removed
    It returns the indexes of the column pairs with value =1 in binary
    matrix, in an MSA that has undergone deletion of several columns;
    the lists gap_list and extreme_list contains the indexes of the
    columns that were deleted, before deletion."""

    logging.info("Retrieve the original coordinates of the matrix.")

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

    logging.info("Retrieving all positions with original coordinates.")

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
    fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
     - %(message)s"""
    logging.basicConfig(filename='mutual_information.log', level=logging.DEBUG,
                        format=fmt)

    msg = 'Runs the computations related to zMIc'
    args_helper = argparse.ArgumentDefaultsHelpFormatter
    argparser = argparse.ArgumentParser(description=msg,
                                        formatter_class=args_helper)
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
    # Prepare the alignment for MIc computations:
    # prune high and low entropy columns
    edited = prune_id_gaps(alignment, args.id)
    gapped_list = []
    if args.g:
        gapped_list = get_all_gaps(edited)
    edited = prune(edited, gapped_list)
    (minlist, maxlist) = get_extreme_columns(edited, args.low, args.high,
                                             args.b)
    edited = prune(edited, minlist+maxlist)
    # compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
    MI_matrix = mutual_info_matrix(edited, args.b)
    ncps_array = NCPS_matrix(edited, args.b)
    MIc_matrix = MI_matrix - ncps_array
    zMIc_matrix = standardise_matrix(MIc_matrix)
    # plot MIc Z-scores and its associated level matrix
    title_zmic = 'zMic of the file {}'.format(args.i)
    plots.plot_heatmap(zMIc_matrix, args.i, title_zmic, args.low)
    tmatrix = get_level_matrix(zMIc_matrix, args.L)
    title_zmic_b = "zMic contacts  with L>{} of the file {}".format(args.L,
                                                                    args.i)
    plots.plot_matrix_binary(tmatrix, args.i, title_zmic_b, args.L)
    # Retrieve CM residue pairs in their original coordinates
    cm_residue_pairs = retrieve_residue_positions(tmatrix, gapped_list,
                                                  minlist + maxlist)
    fd = open(args.o, "w")
    for residue_pair in sorted(cm_residue_pairs):
        fd.write("%d %d\n" % residue_pair)
    fd.close()

    logging.info("Ended program")
    exit(0)
