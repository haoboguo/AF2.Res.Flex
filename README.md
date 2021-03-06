# AlphaFold2 anticipates residue flexibility encoded in protein primary sequences

This repository provides protein structures constructed by AlphaFold2 (AF2), and the codes used in the manuscript "AlphaFold2 models indicate that protein sequence determines both structure and dynamics", which shows that the AF2 models can decipher dynamics information from the protein primary sequences; contributed by Guo et al. (https://www.nature.com/articles/s41598-022-14382-9)

## The AF2 models for 11 proteins (also shown in Table 1; sequences in Appendix)

A: Lanmodulin; 

B: dehalogenase from Deftia acidovorans; 

C: the PAS-A domain protein; 

D: an antifreeze protein (type III) with the X-ray crystallographic structure avalable (PDB ID 1HG7); 

E: a two-domain proten GNE from Homo sapiens (UniProt ID: Q9Y223); 

F: the PAS-A domain containing kinase from Homo sapiens (large protein, the UniProt ID: Q96RG2);

G: the inaZ ice nucleation protein from Pseudomonas syringae (large protein, the UniProt ID P06620); 

H: a heterodimer with the PAS-A domain and the kinase domain; 

I: a homodimer of a MerR-family protein from Mycobacterium tuberculosis (UniProt ID O53384);

J: an intrinsically disordered protein NVJP-1 from Nereis virens; 

K: a randomized protein.


## Python codes for un-pickle the data from AF2 models

plddt.py: get the PLDDT scores for the best AF2 model

pae.py: get the PAE matrix for the best AF2 model

The best model ("ranked_0") that has the highest mean pLDDT score will be pickled and the codes can be modified for unpicling other models.


## R codes (using the GNE protein as examples)

pae.R: the PAE heatmap from the predicted aligned error matrix providied by AF2 (the csv file produced by the python unpicle code will be used as input).

dist.mat2.R: the distance variation matrix estimated and plotted from the MD trajectory (100 ns). A protein pdb ("protein.pdb") and a dcd traject file ("protein.rst1.dcd") will be required as inputs.

pca.R: principal component analysis (PCA) of from MD trajectory (100 ns). The residue cross-correlation analysis is also performed; pymol is required for visualization (see the manuscript). A protein pdb and a dcd trajectory file will be required as input. Here, the PCA will be performed on the backbone atoms (C, O, N, CA); note that the trajectory should contains at least 3N frames (N is the total atom number used in the PCA). Each MD trajectory used in this MS has 10k frames (trajectories have been saved every 10 ps).

The R package "bio3d" is used for trajectory analyses; the "heatmap.2" function in package "gplots" is used to plot the heatmaps. The color scheme is white for high and darkgreen for low PAE/DV scores, respectively, consistent with the schme used in the AF2 database.


## The primary movement (PC1 from the principal component analysis of a 100 ns MD) examples

The two domain protein GNE (system E)
![PC1, sysE](E.pc1.gif)

The homodimer, MerR-family protien (system I)
![PC1, sysI](I.pc1.gif)
