#!/bin/bash
###### USERS INPUT ############################################################

#fluid properties
#Phase 1 is water, phase 2 CO2
#fluid viscosities (m2/s)
Visc1=1e-6
Visc2=8e-06

#fluid densities (kg/m3)
rho1=1000
rho2=1.87

#interfacial tension (N/m)
ift=0.05

#### END OF USER INPUT #######################################################

set -e

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

cp -r 0_org 0
rm 0/Species
sed -i "s/Uinj/0.0/g" 0/U

dimensions=$(cat system/dimensions)
IFS=' ' read -r -a dim_array <<< "$dimensions"
lx1=${dim_array[0]}
lx2=${dim_array[1]}
ly1=${dim_array[2]}
ly2=${dim_array[3]}
lz1=${dim_array[4]}
lz2=${dim_array[5]}

cp system/setFieldsDict1 system/setFieldsDict
sed -i "s/lx1/$lx1/g" system/setFieldsDict
sed -i "s/lx2/$lx2/g" system/setFieldsDict
sed -i "s/ly1/$ly1/g" system/setFieldsDict
sed -i "s/ly2/$ly2/g" system/setFieldsDict
sed -i "s/lz1/$lz1/g" system/setFieldsDict
sed -i "s/lz2/$lz2/g" system/setFieldsDict

echo "setFields"
setFields > setFields0.out

if [ -d "processor0" ]
then
    echo "decomposePar"
    decomposePar -fields > decomposePar0.out

    rm -rf 0
fi
