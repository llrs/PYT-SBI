# cozmic
======
Correlated mutations and distance correlations to predict Aa interactions

Cozmic.py is a stand-alone program which has been fully developed and tested using Python version 3.4.3 and compatible modules.


## Dependencies

The program relays on [biopython](www.biopython.org), [numpy](www.numpy.org), [matplotlib](www.matplotlib.org), and [modeller](https://salilab.org/modeller/).

## Tutorial

The tool can be used in two different ways: either as a full workflow or module by module. Used as a full workflow the user provides a PDB id or file alongside the options of the experiment that must be conducted, then the program returns the main outputs and logging files produced throughout. Used module by module, the user can run the different processes one at a time and also skip those parts which the user is able to solve elsewhere.

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
