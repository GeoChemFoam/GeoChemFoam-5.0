#!/bin/bash

###### USERS INPUT ############################################################

## Define inlet veclocity (m/s)
Ux=0.0001
Uy=0
Uz=0

## Define the kinematic viscocity of the fluid (m^2/s)
##(e.g for water this is 1e-6, for air this would be 1.478e-5)
Visc=1e-06


#### END OF USER INPUT #######################################################

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties

cp -r 0_orig 0
sed -i "s/Ux/$Ux/g" 0/U
sed -i "s/Uy/$Uy/g" 0/U
sed -i "s/Uz/$Uz/g" 0/U

if [ -d "processor0" ]
then
    if [ ! -d "constant/polyMesh" ]
    then
            echo -e "reconstruct parallel mesh"
            reconstructParMesh -constant > reconstructParMesh.out
    fi

    echo -e "DecomposePar"
    decomposePar -fields > decomposeParFlow.out
    rm -rf 0
fi

echo -e "Case initialised"
