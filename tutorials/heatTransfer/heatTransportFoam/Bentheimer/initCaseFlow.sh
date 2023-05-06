#!/bin/bash

###### USERS INPUT ############################################################

#Define flow rate
flowRate=4.2e-9

#fluid properties
Visc=8.2e-8

#Kozeny-Carman constant
kf=1.8e12

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

mkdir -p 0

cp 0_orig/U 0/.
cp 0_orig/p 0/.
sed -i "s/flow_rate/$flowRate/g" 0/U


if [ -d "processor0" ]
then
    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParFlow.out

    rm -rf 0
fi 

echo -e "Case initialised"




