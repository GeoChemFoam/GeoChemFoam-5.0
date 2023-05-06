#!/bin/bash

###### USERS INPUT ############################################################

## Define the diffusion coefficient of the species as a solute (m^2/s)
DiffAB=1e-9
DiffA=1e-9
DiffB=1e-9

#### END OF USER INPUT #######################################################

cp constant/thermoPhysicalProperties1 constant/thermoPhysicalProperties
sed -i "s/DiffAB/$DiffAB/g" constant/thermoPhysicalProperties
sed -i "s/DiffA/$DiffA/g" constant/thermoPhysicalProperties
sed -i "s/DiffB/$DiffB/g" constant/thermoPhysicalProperties

if [ -d "processor0" ]
then
    rm -rf processor*/0/uniform
    
    mkdir 0
    cp 0_orig/A 0/.
    cp 0_orig/B 0/.
    cp 0_orig/AB 0/.

    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParRT.out

    rm -rf 0
else
    rm -rf 0/uniform
    cp 0_orig/A 0/.
    cp 0_orig/B 0/.
    cp 0_orig/AB 0/.
fi

