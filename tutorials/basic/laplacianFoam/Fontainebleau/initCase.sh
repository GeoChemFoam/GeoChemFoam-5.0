#!/bin/bash

###### USERS INPUT ############################################################

## Define your pressure drop at the inlet (Pa)
deltaT=1

## Define diffusion coefficient (m^2/s)
Diff=1e-6

#### END OF USER INPUT #######################################################

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/transportProperties

cp -r 0_orig 0
sed -i "s/deltaT/$deltaT/g" 0/T

if [ -d "processor0" ]
then

    if [ ! -d "constant/polyMesh" ]
    then
            echo -e "reconstruct parallel mesh"
            reconstructParMesh -constant > reconstructParMesh.out
    fi
    echo -e "decomposePar"
    decomposePar -fields > decomposeParT.out
    rm -rf 0
fi

echo -e "Case initialised"

