#!/bin/bash

###### USERS INPUT ############################################################

#Define flow rate
flowrate=5.2075e-11

#fluid viscosities (m2/s) 
#Phase 1 is resident phase, phase 2 is injected phase
Visc1=1e-06
Visc2=7.74e-6

#fluid densities (kg/m3)
rho1=1000
rho2=840

#interfacial tension (N/m)
ift=0.01

#contact angle (degree)
theta=135

#iSpecies molecular diffusion in water (m2/s)r
Diff=1e-9

#Inleti concentration
Ca2plus=8.984e-4
Clminus=1.7968e-3
Hplus=1.051e-7
OHminus=1.069e-7

#Surface site density (kmol/m2)
Surf_aminus=2.0167e-10
Surf_aCaplus=1.7843e-10
Surf_aH=2.0199e-9

#### END OF USER INPUT #######################################################

cp constant/transportProperties1 constant/transportProperties
cp constant/thermoPhysicalProperties1 constant/thermoPhysicalProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties

cp -r 0_org 0

sed -i "s/contactAngle/$theta/g" 0/alpha.water
sed -i "s/flowrate/$flowrate/g" 0/U
sed -i "s/cinlet/$Ca2plus/g" 0/Ca+2
sed -i "s/cinlet/$Clminus/g" 0/Cl-
sed -i "s/cinlet/$Hplus/g" 0/H+
sed -i "s/cinlet/$OHminus/g" 0/OH-

sed -i "s/SurfDen/$Surf_aminus/g" 0/Surf_a-
sed -i "s/SurfDen/$Surf_aCaplus/g" 0/Surf_aCa+
sed -i "s/SurfDen/$Surf_aH/g" 0/Surf_aH


if [ -d "processor0" ]
then
    if [ ! -d "constant/polyMesh" ]
    then
            echo -e "reconstruct parallel mesh"
            reconstructParMesh -constant > reconstructParMesh.out
    fi

    echo -e "Decompose parallel mesh"
    decomposePar -fields > decomposeParRT.out

    rm -rf 0
fi

echo -e "Case initialised"
