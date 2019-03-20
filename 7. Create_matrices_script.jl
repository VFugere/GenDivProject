######################################################################## # Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Science, 2016
# Code in this file by Sen Li and Michael K. Borregaard
#########################################################################

# This code can be run if you have .coords files in your directory
#grid_coordinates("rbcL_noNA", 5/60, 5/60, (0,0))  #apply 4x4 grid cells to the coordinates before calculating GD
#create_master_matrices("rbcL_noNA", "gridcells")

#apply 10-degree wide latitudinal bands to the coordinates
#latband_coordinates("mammals_coi", 10)    # apply 10-degree wide latitudinal bands before calculating GD
#create_master_matrices("mammals_coi", "latbands")

# and so forth if you have data for biomes, equal-area-cells etc
#create_master_matrices("mammals_coi", "anthromes")
#create_master_matrices("mammals_coi", "biomes")
#create_master_matrices("mammals_coi", "equal_area")

include("Create_Master_Matrices_0.5.jl")
create_master_matrices("birds", "coords")
