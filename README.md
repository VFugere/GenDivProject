# GenDivProject

Code for data manipulation and analysis for Millette, Fugère, Debyser et al., "No consistent effects of humans on animal genetic diversity worldwide", _Ecol Lett_

Files numbered 1-7 are bash, R, or Julia scripts used to assemble our dataset, downloading FASTA and .coords files from data repositories, clean them up, align sequences, and compute sequence pairwise comparions. Coords and FASTA files are archived on Dryad (see link provided in manuscript).

Files 8-12 calculate the mean pairwise geographic distance among sequences (D), divide files into populations, compute population genetic diversity from mean pairwise sequence dissimilarity (pi), load and format HYDE 3.2 data files, and assign land use and human density values to all sequences and populations. The product of these scripts are the two formatted datasets used in all analyses (sequence_metadata.RData & DF_Master.RData; both of these are available in the /Data folder).

Files 13-19 use the two formatted datasets to fit GAMMs and to produce the figures of the manuscript.
