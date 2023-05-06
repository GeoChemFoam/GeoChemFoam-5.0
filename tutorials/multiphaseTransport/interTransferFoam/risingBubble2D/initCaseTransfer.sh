#!/bin/bash
###### USERS INPUT ############################################################

#fluid properties
#Phase 1 is resident fluide, phase 2 is the bubble
#fluid viscosities (m2/s)
Visc1=2e-5
Visc2=1.5e-05

#fluid densities (kg/m3)
rho1=1200
rho2=1.2

#interfacial tension (N/m)
ift=0.065

#Species properties
#Diffusion in each phase (m2/s)
Diff1=2e-5
Diff2=0 #the gas is pure
#Henry's constant cphase2=H*cphase1
H=10
#molecular weight (kg/kmol)
Mw=44

#### END OF USER INPUT #######################################################

cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

cp constant/thermoPhysicalProperties0 constant/thermoPhysicalProperties
sed -i "s/Diff_1/$Diff1/g" constant/thermoPhysicalProperties
sed -i "s/Diff_2/$Diff2/g" constant/thermoPhysicalProperties
sed -i "s/H_12/$H/g" constant/thermoPhysicalProperties
sed -i "s/M_w/$Mw/g" constant/thermoPhysicalProperties

if [ -d "processor0" ]
then
    reconstructPar -withZero > reconstructParTransfer.out
fi

SpeciesMol=$(echo "scale=13; $rho2/$Mw" | bc)
rm -f 0/Species
rm -f tmp tmp2
awk 'NR < 24 { print } ' 0/alpha.phase1 > 0/Species
sed -i "s/0 0 0 0 0 0 0/0 -3 0 0 1 0 0/g" 0/Species
sed -i "s/alpha.phase1/Species/g" 0/Species
awk 'NR > 23 { print } ' 0/alpha.phase1 > tmp
head -n -21 tmp > tmp2
awk '{ print ('$SpeciesMol')*(1-$l)}' tmp2 >> 0/Species
tail -21 0/alpha.phase1 >> 0/Species
rm -f tmp tmp2 

cp constant/g1 constant/g

if [ -d "processor0" ]
then

    echo -e "decompose"
    decomposePar -fields > decomposeParT.out

    rm -rf 0
fi




