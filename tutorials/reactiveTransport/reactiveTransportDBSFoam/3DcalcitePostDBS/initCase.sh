#!/bin/bash

###### USERS INPUT ############################################################

#Define flow rate
flowRate=3.5e-10

#fluid properties
Visc=2.61e-6
Diff=5e-9

#Reaction constants
kreac=8.9125e-4 
scoeff=2
rhos=2710
Mws=100
cinlet=0.0126


#Kozeny-Carman constant
kf=1.8e12 

#### END OF USER INPUT #######################################################

cp constant/dynamicMeshDict0 constant/dynamicMeshDict
cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution

echo -e "set flow and transport properties"
cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/rho_s/$rhos/g" constant/transportProperties
sed -i "s/Mw_s/$Mws/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

cp constant/thermoPhysicalProperties0 constant/thermoPhysicalProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties
sed -i "s/s_coeff/$scoeff/g" constant/thermoPhysicalProperties
sed -i "s/k_reac/$kreac/g" constant/thermoPhysicalProperties

if [ -d "processor0" ]
then
    mkdir 0
    cp 0_orig/U 0/.
    cp 0_orig/C 0/.
    cp 0_orig/p 0/.

    sed -i "s/flow_rate/$flowRate/g" 0/U
    sed -i "s/c_inlet/$cinlet/g" 0/C

    echo "decompose parallel mesh"

    decomposePar -fields > decomposeParFields0.out
    rm -rf 0
else
    cp 0_orig/U 0/.
    cp 0_orig/C 0/.
    cp 0_orig/p 0/.

    sed -i "s/flow_rate/$flowRate/g" 0/U
    sed -i "s/c_inlet/$cinlet/g" 0/C
fi

echo "Case initialised"
