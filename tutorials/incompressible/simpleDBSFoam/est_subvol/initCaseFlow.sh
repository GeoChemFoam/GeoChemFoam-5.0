#!/bin/bash

###### USERS INPUT ############################################################

## Define your pressure drop at the inlet (Pa)
PGRAD=100000

## Define the kinematic viscocity of the fluid (m^2/s)
##(e.g for water this is 1e-6, for air this would be 1.478e-5)
Visc=1e-6

#### END OF USER INPUT #######################################################

cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties

cp constant/fvOptionsRun constant/fvOptions
sed -i "s/PGRAD/$PGRAD/g" constant/fvOptions

if [ -d "processor0" ]
then
    cp -r 0_orig 0
    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParFlow.out

    rm -rf 0
else
    cp 0_orig/* 0/.
fi

echo "Case initialised"
