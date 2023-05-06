#!/bin/bash

###### USERS INPUT ############################################################

#Define flow rate
PDROP=0.00005

#fluid properties
Visc=1e-6

#Kozeny-Carman constant
kf=1.8e8

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

cp system/fvSolution1 system/fvSolution

if [ -d "processor0" ]
then
    mkdir 0
    cp 0_orig/U 0/.
    cp 0_orig/p 0/.
    sed -i "s/PDROP/$PDROP/g" 0/p

    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParFlow.out

    rm -rf 0
else
    cp 0_orig/U 0/.
    cp 0_orig/p 0/.
    sed -i "s/PDROP/$PDROP/g" 0/p
fi

echo -e "Case initialised"




