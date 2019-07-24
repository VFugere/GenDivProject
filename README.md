# GenDivProject

Code for data manipulation and analysis for Gonzalez Lab project on human impacts on animal genetic diversity.

Files numbered 1-12 are bash, R, or Julia scripts used to assemble our dataset, downloading FASTA and .coords files from data repositories, clean them up, and use them to calculate population genetic diversity (mean pairwise nucleotide diversity) and spatial distance (mean pairwise geographic distance among sequences). Files 10-11 format HYDE 3.2 data files and assign land use and human density values to all sequences. File 12 merges genetic distance, geographic distance, and land use datasets into a master dataset used in analyses.

The two datasets necessary to run data visualization and analysis scripts (Files 13-18) are available in the /Data folder. 
