# cozmic

Correlated mutations and distance correlations to predict Aa interactions

Cozmic.py is a stand-alone program which has been fully developed and tested using Python version 3.4.3 and compatible modules.


## Dependencies

The program relays on [biopython](www.biopython.org), [numpy](www.numpy.org), [matplotlib](www.matplotlib.org), and [modeller](https://salilab.org/modeller/).

## Tutorial

The tool can be used in two different ways: either as a full workflow or module by module. Used as a full workflow the user provides a PDB id or file alongside the options of the experiment that must be conducted, then the program returns the main outputs and logging files produced throughout. Used module by module, the user can run the different processes one at a time and also skip those parts which the user is able to solve elsewhere.

### cozmic

This is the main module. The purpose of the program is to study protein contact and contact prediction carried out with zMIc, an information-theoretic method that ultimately resorts to an MSA aligned to the target structure.

**Input and Options**: In all cases the input shall be a PDB structure. The user mustdecide which type of analysis wants to perform: analysis of known PDB structure or derive a model from the input. One possible motivation to study a model created from a PDB might be to carry out comparative analysis to test the robustness of contact predictions with proteins that have the same underlying sequence, but perturbed structures. The types of analysis are under the subcommand real and model, respectively. An optional argument set to tune the verbosity of the logger has also been added.

Under the subcomand **real** we have the following arguments and options available:
+ Required arguments:
  + input: either the id of a PDB or the path to a PDB file
+ Optional arguments:
  + -a: Atom to calculate distance between residues of the structure, can be either CA, CB or min, see module contact map
  + -CA: Threshold for Cα distance.
  + -CB: Threshold for Cβ distance.
  + -min: Threshold for minimum distance between residues.
  + -blast: Type of BLAST the user wants to perform for retrieving related sequences of the PDB. It can be: blastp, tblastn, blastx, tblastx.
  + -db: Database against which the BLAST is done. It can be: pdb, swissprot, nr, tsa nr, env nr.
  + -s: Maximum number of hits resulting from the BLAST.
  + -f: Filter the output of BLAST by genus or not (default yes).
  + -g: Delete columns of the MSA with at least one gap for the calculus of MI.
  + -b: Base of the logarithm.
  + -low: Minimum entropy for the computation of MI.
  + -high: Maximum entropy for the computation of MI.
  + -m: Program to produce the MSA options are [muscle](http://www.ebi.ac.uk/Tools/msa/muscle/), [t-coffe](http://tcoffee.crg.cat/) or [clustalw](http://www.genome.jp/tools/clustalw/).

Under the **model** subcommand 
+ Required arguments:
  + seq: Sequence to modelize
  + models: Sequence which is the template
+ Optional arguments
  + -pir: Alignment in pir format required
  + -fasta: Alignment in fasta format
  
**Output**: The program produces several plots and logging files: distance map plot, contact map plot, plot zMIc values, plot with contact predictions based on zMIc, plot with number of predicted contacts and precision of prediction against the cut-off level employed to select candidate residue pairs, file with BLAST outcome after non-redundancy filtering, file with MSA outcomes, file enclosing the list of residue pair coordinates which have zMIc ≥ 2.

#### Examples
##### Example 1
In this example we will run cozmic to analyse the contact predictions for the PDB structure 1cd8. We want the program to count distances between C β atoms. We also want it to build the MSA based on at most 1000 BLASTP hits filtered by genus to
prevent redundancy artifacts. We commit to use Muscle as MSA method of choice.

`> python3 cozmic.py real 1cd8 -a CB -s 1000 -m muscle`

Depending on the remote availablity of NCBI’s QBl A st service and the computational load the MSA has to deal with, the program may last for about 1-2 minutes to provide all the results. At the end of the execution, the following new files would have been generated:

Graphics files: distance map, contact map (see Figure 6), zMIc heatmap, zMIc predicted contacts and zMIc prediction precision diagram, enclosing the curves of number of predicted hits and precision versus zMIc cut-off, respectively.
```
contact_map_1cd8_c_CB.png
heatmap_1cd8_d_CB.png
heatmap_1cd8_z_0.3.png
contact_map_1cd8_p_CB.png
cutoff_zMIc_p_1cd8.png
```

Text files: query passed to BLAST, result from doing BLAST after non-redundancy filtering, MSA obtained and log including the residue pairs that have zMIc value ≥ 2.
```
cozmic_blast_query.fa
cozmic_blast.out
aligned.fa
cozmic.log
```

##### Example 2
In this example we will run cozmic to analyse the contact predictions for the PDB structure 1akj. We want the program to count distances between C α atoms and we want it to consider contact as being at a distance ≤ 15 Angstroms. We also want it to build the MSA based on at most 2000 BLASTP hits filtered by genus to prevent redundancy artifacts.
We first commit to use Muscle as an MSA method of choice.

`> python3 cozmic.py real -s 2000 -a CA -CA 15 1akj -m muscle`

Recall that when executing this command, the default option set for the BLAST database is nr (non-redundant database) and the non-redundancy filtering by genus is set by default.

Depending on the remote availablity of NCBI’s QBl A st service and the computational load the MSA has to deal with, the program may last for about several minutes to provide all the results. At the end of the execution, the following new files would have been generated:

Graphics files: distance map, contact map, zMIc heatmap, zMIc predicted contacts and zMIc prediction precision diagram.
```
contact_map_1akj_c_CA.png
heatmap_1akj_d_CA.png
heatmap_1akj_z_0.3.png
contact_map_1akj_p_CA.png
cutoff_zMIc_1akj.png
```

Text files: query passed to BLAST, result from doing BLAST after non-redundancy filtering, MSA obtained and log including the residue pairs that have zMIc value ≥ 2.
```
cozmic_blast_query.fa
cozmic_blast.out
aligned.fa
cozmic.log
```

### contact_map
The module contact map provides the functions needed to produce a contact map based on the protein’s structure information given in a PDB file. It uses the standard Aa’s to calculate the pairwise distances between residues. All the residues of any chain of all models in the PDB are used. Both heteroatoms and ligands are omitted from distance
computations.
Remark: Although PDB files usually include the positional information of regions ofproteins, they may not account for all the residues of the protein due to the experimental difficulties to ascertain certain atomic positions, particularly for highly flexible regions.

**Input**: The input is the first argument and should be the structure of a protein in PDB format.

**Options**: The module can work under three conditions set by command line after the **-a** option:
 The default option measures the distance between residues as the minimum distance
between non-hydrogen atoms of the respective residues. Two residues are labelled
as contact if their distance is below 6 Angstroms.
+ If the option -a CA is set, then only those residues whose respective Cα are closer
than 16 Angstroms are taken into account.
+ If the option **-a** CB is set, then only those residues whose respective Cβ are closer are taken into account. **Note**: The amino acid glycine does not have any Cβ, In this case this option the distance is taken from its Cα.
There is also the possibility to set the threshold for each distance defintion:
+ With the option **-CA** the threshold value for Cα distances can be specified.
+ With the option **-CB** the threshold value for Cβ distances can be specified.
+ With the option **-min** the threshold value for minimal distances can be specified.

*Remark*: Those residue pairs with a positional difference in the primary structure d ≤ 3
are excluded from being considered contact.

**Output**: The module produces two plots:
+ A heatmap representing the distances for each pair of residues in the structure.
This plot is stored in a `heatmap_{pdbname}_{option}.png` file.
+ A plot representing the the pairs of residues labelled as contact. This plot is stored
as contact `map_{pdbname}_{option}.png`. The balck spots represent contacts
between residues, whence a putative interaction.

*Remark*: Both plots arise from the known structure of the protein encoded as a PDB;
consequently, the indexes used to label residues refer to the PDB structure, not to the
protein sequence.

#### Examples

`> python3 contact_map.py pdb1cd8.ent`
Will produce the following files 
```
distance_map_pdb1cd8_min.png
contact_map_pdb1cd8_min.png
contact_map.log
```
Or you can specify that you want to calculate the Ca distance at 10 Angstroms with:

`> python3 contact_map.py pdb1cd8.ent -a CA -CA 10`
Will produce the following files 
```
distance_map_pdb1cd8_CA.png
contact_map_pdb1cd8_CA.png
contact_map.log
```

