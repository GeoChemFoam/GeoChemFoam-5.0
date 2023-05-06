#!/bin/bash

###### USERS INPUT ############################################################

## Define your pressure drop at the inlet (Pa)
deltaP=10


## Define the kinematic viscocity of the fluid (m^2/s)
##(e.g for water this is 1e-6, for air this would be 1.5e-5)
Visc=1.5e-5

#### END OF USER INPUT #######################################################

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties

mkdir 0
cp -r 0_orig/U 0/.
cp -r 0_orig/p 0/.
sed -i "s/deltaP/$deltaP/g" 0/p

if [ -d "processor0" ]
then

    if [ ! -d "constant/polyMesh" ]
    then
            echo -e "reconstruct parallel mesh"
            reconstructParMesh -constant > reconstructParMesh.out
    fi
    echo -e "decomposePar"
    decomposePar -fields > decomposeParFlow.out
    rm -rf 0
fi

echo -e "Case initialised"

