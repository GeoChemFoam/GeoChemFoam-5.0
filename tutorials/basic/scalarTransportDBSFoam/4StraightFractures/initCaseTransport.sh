#!/bin/bash

###### USERS INPUT ############################################################

#Diffusion coefficient (m^2/s)
Diff=1e-9

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties2 constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/transportProperties

if [ -d "processor0" ]
then
    mkdir 0
    cp -r 0_orig/T 0/.
    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParT.out
    rm -rf 0
else
    cp -r 0_orig/T 0/.
fi

echo -e "Case initialised"


