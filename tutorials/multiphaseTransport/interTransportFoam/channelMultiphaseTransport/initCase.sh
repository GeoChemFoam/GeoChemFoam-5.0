#!/bin/bash
###### USERS INPUT ############################################################

## Phase and species initialised as a square
x1=0
x2=4e-5
y1=0
y2=5e-5

#fluid properties
#Phase 1 is injected phase, phase 2 is resident phase
#fluid viscosities (m2/s)
Visc1=18e-6
Visc2=1.52e-06

#fluid densities (kg/m3)
rho1=1
rho2=789

#interfacial tension (N/m)
ift=0.02

#Species properties
#Diffusion in each phase (m2/s)
Diff1=2e-7
Diff2=2e-7
#Henry's constant cphase2=H*cphase1
H=0.1

#Injection velocity
Ux=0.4

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
cp -r 0_orig 0
sed -i "s/Ux/$Ux/g" 0/U


cp system/setFieldsDict1 system/setFieldsDict
sed -i "s/x1/$x1/g" system/setFieldsDict
sed -i "s/x2/$x2/g" system/setFieldsDict
sed -i "s/y1/$y1/g" system/setFieldsDict
sed -i "s/y2/$y2/g" system/setFieldsDict

echo "set fields"
setFields > setFieldsRun.out

if [ -d "processor0" ]
then
    echo "decomposePar fields"
    decomposePar -fields > decomposeParRun.out

    rm -rf 0
fi

echo -e "Case initialised"


