# AF2.Res.Flex
AlphaFold2 Predicts Residue Flexibility

This repository contains codes for analysis of the AlphaFold2 (AF2) predicted aligned error (PAE) maps, and the distance deviation (DV) maps from molecular dynamics (MD) trajectories. We show that AF2 accurately predict the residue flexibilities captured by MD simulations. Code of the principal component analysis (PCA) of MD trajectories is also provided, together with the primary movements (PC1) of selected systems.

The R packages "bio3d" is used to analyze the MD trajectories (dcd format) for Distance Variation (DV) calculations and PCA analysis; the "ggplot2", "gplots" are used for the plots.

The PAE maps are plotted from the pkl files from the AF2 modeling. The unpickle script and R script for the heatmap are provided. The DV maps are calculated using the alpha-carbon atoms of all residues. The mass-weighted PCA analysis are performed to all backbone heavy atoms.
