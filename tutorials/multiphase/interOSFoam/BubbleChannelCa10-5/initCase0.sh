#!/bin/bash
###### USERS INPUT ############################################################

#Bubble initially defined as rectangle from x1 to x2
x1=0.00002
x2=0.00022

#fluid viscosities (m2/s)
#Phase 1 is resident phase, phase 2 is injected phase
Visc1=1.52e-06
Visc2=1.8e-5

#fluid densities (kg/m3)
rho1=789
rho2=1

#interfacial tension (N/m)
ift=0.02

#contact angle (degree)
theta=0

#### END OF USER INPUT #######################################################

cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

cp -r 0_orig 0

sed -i "s/contactAngle/$theta/g" 0/alpha.phase1
sed -i "s/Ux/0.0/g" 0/U

cp system/setFieldsDict1 system/setFieldsDict
sed -i "s/x1/$x1/g" system/setFieldsDict
sed -i "s/x2/$x2/g" system/setFieldsDict

echo "set initial bubble as rectangle"
setFields > setFields0.out

if [ -d "processor0" ]
then
    echo "decomposePar"
    decomposePar -fields > decomposePar0.out

    rm -rf 0
fi

echo -e "Case initialised"
