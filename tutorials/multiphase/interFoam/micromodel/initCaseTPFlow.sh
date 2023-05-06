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

#### END OF USER INPUT #######################################################

cp constant/transportPropertiesTP constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

rm -rf 0 [1-9]*
rm -rf processor*/0
rm -rf processor*/[1-9]*

cp -r 0_org 0

sed -i "s/contactAngle/$theta/g" 0/alpha.water
sed -i "s/flowrate/$flowrate/g" 0/U

echo -e "set alpha field"
setFields > setFieldsTP.out

if [ -d "processor0" ]
then
    decomposePar -fields  > decomposeParTP.out
    rm -rf 0
fi


