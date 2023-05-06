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

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties

cp constant/thermoPhysicalProperties0 constant/thermoPhysicalProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties


rm -rf 0_org
cp -r 0_gold 0_org
sed -i "s/flow_rate/$flowRate/g" 0_org/U
sed -i "s/k_reac/$kreac/g" 0_org/C
sed -i "s/s_coeff/$scoeff/g" 0_org/C
sed -i "s/c_inlet/$cinlet/g" 0_org/C
sed -i "s/k_reac/$kreac/g" 0_org/pointMotionU
sed -i "s/Mw_s/$Mws/g" 0_org/pointMotionU
sed -i "s/rho_s/$rhos/g" 0_org/pointMotionU

cp -r 0_org 0


cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution

if [ -d "processor0" ]
then
    if [ ! -d "constant/polyMesh" ]
    then
            echo -e "reconstruct parallel mesh"
            reconstructParMesh -constant > reconstructParMesh.out
    fi

    echo "decompose parallel"
    decomposePar -fields > decomposePar0.out
    for i in processor*;do cp -r $i/constant/polyMesh $i/0/.; done

    rm -rf 0
else
    cp -r constant/polyMesh 0/.
fi

echo -e "Case initialised."
