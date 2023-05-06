# Automation of Immune Peptide Generation

## Overview

The `bioinf.ipynb` file contains a script divided into blocks that allows the design of potentially immunogenic peptides from a protein structure. The script loads the protein structure from a file, identifies the cleavage site, and finds the available region in the protein. Then, candidate peptides of different lengths in this region are generated, their immunogenicity is calculated, and already known antigenic motifs are filtered. Next, peptide sequences are optimized based on their immunogenicity, and the best peptides are selected for display.

## Parameters

The script allows for the variation of the following parameters:
- `file_path` - path to the protein structure file.
- `site_sequence` - amino acid sequence to be cleaved for peptide generation (taken from the `art.pdf` article).
- `sasa_threshold` - threshold value for calculating available regions in the protein structure.
- `peptide_length_range` - range of possible peptide lengths.
- `known_antigenic_motifs` - known antigenic motifs to be excluded from peptide generation.
- `top_n` - number of best peptides to be selected for optimization.

## In silico Molecule Analysis

For further analysis of linear peptide sequences, I propose using AlphaFold2 to convert them to the pdb format. Then, through obabel, pdb can be converted to pdbqt and the resulting files can be loaded into the `docking/ligands` folder. (this can be done as a script)

The protein structure in pdbqt format should be placed in the `receptors` folder.
The peptide structures in pdbqt format should be placed in the `ligands` folder.
Then you can run the Vina script, which will create configuration files for AutodockVina and run molecular docking. Manual adjustment may be required.

The obtained results can be analyzed, and additional research on molecular dynamics can be conducted.
