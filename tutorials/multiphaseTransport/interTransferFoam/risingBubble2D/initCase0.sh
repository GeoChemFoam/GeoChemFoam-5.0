#!/bin/bash
###### USERS INPUT ############################################################

## Phase and species initialised as a square
x1=0
x2=0.002
y1=0.002
y2=0.0042

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

#### END OF USER INPUT #######################################################

cp constant/transportProperties0 constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

cp -r 0_orig 0

cp system/setFieldsDict1 system/setFieldsDict
sed -i "s/x1/$x1/g" system/setFieldsDict
sed -i "s/x2/$x2/g" system/setFieldsDict
sed -i "s/y1/$y1/g" system/setFieldsDict
sed -i "s/y2/$y2/g" system/setFieldsDict

echo "initialise bubble as square"
setFields > setFields0.out

cp constant/g0 constant/g

if [ -d "processor0" ]
then 
    echo -e "decompose"
    decomposePar -fields > decomposePar0.out

    rm -rf 0
fi

echo -e "Case initialised"
