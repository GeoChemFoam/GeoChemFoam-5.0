#!/bin/bash

###### USERS INPUT ############################################################

#Define flow rate
flowrate=3.5e-8

#fluid viscosities (m2/s) 
#Phase 1 is resident phase, phase 2 is injected phase
Visc1=1e-06
Visc2=1.65e-5

#fluid densities (kg/m3)
rho1=1000
rho2=864

#interfacial tension (N/m)
ift=0.03

#contact angle (degree)
theta=45

#Phase-field interface thickness parameters (at least 1.5x larger than mesh resolution)
eps=1e-4

#Phase-field mobility parameter (set up so that local Peclet number ~10) 
mob=2e-7

#### END OF USER INPUT #######################################################

cp constant/transportPropertiesTP constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties
sed -i "s/thick/$eps/g" constant/transportProperties
sed -i "s/mobi/$mob/g" constant/transportProperties

rm -rf 0 [1-9]*
rm -rf processor*/0
rm -rf processor*/[1-9]*

cp -r 0_org 0

dimensions=$(cat system/dimensions)
IFS=' ' read -r -a dim_array <<< "$dimensions"
x_dim=${dim_array[0]}
y_dim=${dim_array[1]}
z_dim=${dim_array[2]}
n_x=${dim_array[3]}
n_y=${dim_array[4]}
n_z=${dim_array[5]}
direction=${dim_array[6]}
l1=${dim_array[7]}

python system/createCEquil.py --xMin 0 --xMax $x_dim --yMin 0 --yMax $y_dim --zMin 0 --zMax $z_dim --direction $direction --linter $l1 --nX $n_x --nY $n_y --nZ $n_z --eps $eps

sed -i "s/contactAngle/$theta/g" 0/C

sed -i "s/flowrate/$flowrate/g" 0/U

if [ -d "processor0" ]
then
    echo -e "Decompose parallel mesh"
    decomposePar -fields > decomposeParTP.out
fi
