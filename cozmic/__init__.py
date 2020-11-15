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
import logging
import os

# Biopython modules
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Data import SCOPData
from Bio import AlignIO

# Non-standard modules
import contact_map as cm
import mutual_information as mut
import msa_caller as msa
import blast as blst
import plots
# import modeller_caller as mc
# from modeller import environ

# Entrez inputs
Entrez.email = "ferran.muinos@gmail.com"
Entrez.tool = "cozmic.py"

# This program may well be called cozmic.py

if __name__ == '__main__':
    msg = 'Runs the main workflow'
    default_help = argparse.ArgumentDefaultsHelpFormatter
    argparser = argparse.ArgumentParser(description=msg,
                                        formatter_class=default_help)
    msg_subparser = 'Choose between modelize a pdb or from an existing file'
    subparsers = argparser.add_subparsers(help=msg_subparser)
    argparser.add_argument("-d",
                           help="""Decide the degree of the logging. This
                           option can be used whatever other options are.""",
                           action="count",
                           default=1)

    real = subparsers.add_parser('real', help='Make a structure of a real pdb',
                                 formatter_class=default_help)
    model = subparsers.add_parser('model', help='Modelize a pdb',
                                  formatter_class=default_help)
    model.add_argument("seq", help="Name of the sequence to analyse")
    model.add_argument("models", help="Models of the file")
    model.add_argument("-pir", help="Name of the file in pir format")
    model.add_argument("-fasta",
                       help="File with sequences in fasta format")
    real.add_argument("input", help="Target sequence id or filename")
    # argument options for cm module functions
    real.add_argument("-a", help="""Aa atom to calculate distance with:
                      CA (alpha carbon), CB (beta carbon), min (minimum
                      distance between non-hydrogen atom pairs from each
                      residue).""",
                      choices=["CA", "CB", "min"],
                      default="min")
    real.add_argument("-CA",
                      help="""Set the threshold distance between Carbon alpha
                      atoms.""",
                      type=int,
                      default=15)
    real.add_argument("-CB",
                      help="""Set the threshold distance between Carbon beta
                      atoms.""",
                      type=int,
                      default=12)
    real.add_argument("-min",
                      help="""Set the minimal threshold distance between
                      atoms.""",
                      type=int,
                      default=6)
    # argument options for blst module functions
    real.add_argument("-blast", help="""Type of BLAST search to be
                      performed. Notice that blastp, blast and tblastn
                      require a protein sequence, whereas blastx and
                      tblastx require a nucleotide sequence.""",
                      choices=["blastp", "tblastn", "blastx", "tblastx"],
                      default="blastp")
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
    real.add_argument("-f", help="""If present, then don't filter the BLAST
                      output by genus for attaining non-redundancy;
                      otherwise filter by genus.""",
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

    fmt = """%(asctime)s - %(filename)s - %(funcName)s - %(levelname)s
      - %(message)s"""
    logging.basicConfig(filename='cozmic.log', level=int(100/(args.d*10)),
                        format=fmt)

    print("You are currently running the program with: ", args)
    try:
        args.input
    except AttributeError:
        from modeller import environ
        import modeller_caller as mc
        env = environ()  # Some variables needed for the modeller
        modeler = mc.modeller_caller(env)
        # Convert the fasta alignment in pir format
        if not args.fasta and not args.pir:
            raise argparser.error("Required a fasta or a pir alignment")
        elif args.fasta:
            modeler.convert_ali(args.fasta, args.pir)
        modeler.modelize(args.pir, args.seq, args.models)
    else:
        # Retrieve the PDB structure, filter and get sequence
        logging.captureWarnings(True)
        parser = PDBParser(PERMISSIVE=1)
        if os.path.isfile(args.input):
            pdbpath = args.input
        else:
            pdbl = PDBList()
            try:
                pdbpath = plots.pdb_download(args.input, os.getcwd())
            except:
                raise FileExistsError("Make sure your query format is correct")
        structure = parser.get_structure("cozmic_pdb_query", pdbpath)
        residues = cm.filter_residues(structure)
        s = ""
        for residue in residues:
            s += SCOPData.protein_letters_3to1.get(residue.get_resname(), 'X')
        seq = Seq(s, generic_protein)
#         sys.stderr.write("Protein sequence:%s\n" % seq)
        # Compute distances and contact between residues
        dist_matrix = cm.calc_dist_matrix(residues, args.a)
        cont_matrix = cm.contact_map(dist_matrix, args.a)
        # Set up the files for calling BLAST
        blast_query_name = "cozmic_blast_query.fa"
        blast_out_name = "cozmic_blast.out"

        fd = open(blast_query_name, "w")
        fd.write(">%s\n%s" % ("cozmic_blast_query", seq))
        fd.close()
        # Enter the query as first element in the output file
        file_out = open(blast_out_name, "w")
        record = SeqIO.read("cozmic_blast_query.fa", format="fasta")
        SeqIO.write(record, file_out, "fasta")
        # Call run_BLAST and write the output in the output file
        # So far it is not working
        blast_result = blst.run_BLAST(blast_query_name, args.blast,
                                      args.db, args.s)
        ides = blst.analyze_blast_result(blast_result, args.f)
        ids = list(blst.filter_ids(ides, "gi"))
        if len(ids) <= 5:
            msg_err = "Seems that a bad output on blast generated no hits."
            raise ValueError(msg_err)
        print(len(ids), ids)
        SeqIO.write(blst.retrive_sequence(ids), file_out, "fasta")
        file_out.close()
        # MSA: align the query to its homologs with the method of choice
        msa.call_msa_method(args.m, blast_out_name, "aligned.fa", "fasta")
        alignment = AlignIO.read("aligned.fa", "fasta")
        # Prepare the alignment for MIc computations:
        # prune high and low entropy columns
        edited = mut.prune_id_gaps(alignment, "cozmic_blast_query")
        gapped_list = []
        if args.g:
            gapped_list = mut.get_all_gaps(edited)
        edited = mut.prune(edited, gapped_list)
        (minlist, maxlist) = mut.get_extreme_columns(edited, args.low,
                                                     args.high, args.b)
        mm = minlist + maxlist
        edited = mut.prune(edited, mm)

        # compute MI, NCPS, MIc, Z-score MIc + its associated level matrix
        MI_matrix = mut.mutual_info_matrix(edited, args.b)
        ncps_array = mut.NCPS_matrix(edited, args.b)
        MIc_matrix = MI_matrix - ncps_array
        zMIc_matrix = mut.standardise_matrix(MIc_matrix)

        # plot distance, contact, MIc Z-scores and its associated level matrix
        title_dist = 'Distance of the file {}'.format(args.input)
        plots.plot_heatmap(dist_matrix, args.input+"_d", title_dist,
                           args.a, "Angstroms")

        # Plot contacts based on distance
        title_binary = 'Contacts of the file {}'.format(args.input)
        plots.plot_matrix_binary(cont_matrix, args.input+"_c",
                                 title_binary, args.a)

        # Plots zMIc
        title_zmic = 'zMic values of the file {}'.format(args.input)
        plots.plot_heatmap(zMIc_matrix, args.input+"_z",
                           title_zmic, args.low, "zMIc")
        tmatrix = mut.get_level_matrix(zMIc_matrix, 2)

        # Calculates the default predicted contacts and store them in the file

        std_cont = mut.retrieve_residue_positions(tmatrix, gapped_list, mm)
        name_out = "predicted_contacts_{}.out"
        with open(name_out.format(args.input), "w") as out_f:
            out_f.write(repr(std_cont))

        # Plots the default contacts
        title_zmic_b = "Predicted contacts with zMic > 2 of the file"
        plots.plot_matrix_binary(tmatrix, args.input+"_p",
                                 title_zmic_b.format(args.input), args.a)

        # plot level-precision analysis and CM-distance analysis
        cutoff_list, hit_list, precision_list = plots.precision_analysis(
        zMIc_matrix, cont_matrix, gapped_list, mm, 0.0, 3.0, 60)
        plots.plot_twin_curves(cutoff_list, hit_list, precision_list,
                               args.input)

        logging.captureWarnings(False)
