#!/bin/bash

###### USERS INPUT ############################################################

# Define cropping parameters
x_min=0
x_max=200
y_min=0
y_max=150
z_min=0
z_max=256

# number of cells of initial mesh
n_x=100
n_y=75
n_z=128

# Level of refinement - mesh is refined by n_level at the pore surfaces
n_level=1

# define resolution (m)
res=0.0000053

#### END OF USER INPUT #######################################################

# ###### DO NOT MAKE CHANGES FROM HERE ###################################

python writeSpeciesHdf5.py --x_min=$x_min --x_max=$x_max --y_min=$y_min --y_max=$y_max --z_min=$z_min --z_max=$z_max --res=$res --n_x=$n_x --n_y=$n_y --n_z=$n_z --n_level=$n_level 


