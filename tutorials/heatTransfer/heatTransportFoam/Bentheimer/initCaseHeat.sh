#!/bin/bash

###### USERS INPUT ############################################################

#fluid properties
kappa_f=0.04 
kappa_s=3.3

rho_f=300
rho_s=2600

gamma_f=3000
gamma_s=700

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties2 constant/transportProperties
sed -i "s/kappa_s/$kappa_s/g" constant/transportProperties
sed -i "s/kappa_f/$kappa_f/g" constant/transportProperties
sed -i "s/rho_f/$rho_f/g" constant/transportProperties
sed -i "s/rho_s/$rho_s/g" constant/transportProperties
sed -i "s/gamma_f/$gamma_f/g" constant/transportProperties
sed -i "s/gamma_s/$gamma_s/g" constant/transportProperties


mkdir -p 0
cp 0_orig/T 0/.


if [ -d "processor0" ]
then
    # Decompose
    echo -e "DecomposePar"
    decomposePar -fields > decomposeParH.out

    rm -rf 0
fi



