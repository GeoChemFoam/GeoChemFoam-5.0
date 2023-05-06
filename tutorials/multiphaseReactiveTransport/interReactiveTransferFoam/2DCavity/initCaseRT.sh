#!/bin/bash
###### USERS INPUT ############################################################

#fluid viscosities (m2/s)
Visc1=1e-6
Visc2=8e-06

#fluid densities (kg/m3)
rho1=1000
rho2=1.87

#interfacial tension (N/m)
ift=0.05

### Species name (the first species needs to be CO2)
Species=('CO2' 'H+' 'OH-' 'CO3-2' 'HCO3-')

### Diffusion coefficient in water###
### original values multiplied by 100 to reduce computational time
Diff=('1.6e-7' '9.83e-7' '5.27e-7' '0.955e-7' '5.27e-7')

### Henry's constant for CO2 (only CO2 can go in the gas phase)
He=1.25

### Species molecular weight
Mw=('44' '1' '17' '60' '61')

### Water concentration (kmol/m3)
Cmol=('0' '9.97413e-8' '1.00954e-7' '0' '0')

### Injection velocity
Uinj=3e-3

#### END OF USER INPUT #######################################################

set -e

echo $Species
if [ $Species != 'CO2' ]
then
    echo "ERROR: first species should be CO2"
    exit
fi

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc1/$Visc1/g" constant/transportProperties
sed -i "s/Visc2/$Visc2/g" constant/transportProperties
sed -i "s/rho1/$rho1/g" constant/transportProperties
sed -i "s/rho2/$rho2/g" constant/transportProperties
sed -i "s/ift/$ift/g" constant/transportProperties

python system/createThermoPhysicalProperties.py --Species ${Species[@]} --Diff ${Diff[@]} --He $He --Mw ${Mw[@]}

if [ -d "processor0" ]
then
    reconstructPar -withZero > reconstructParRT.out
fi


index=0
for i in ${Species[@]}
do
    cw=${Cmol[$index]}
    if [ $i == 'CO2' ]
    then
        cg=$(echo "scale=13; $rho2/${Mw[0]}" | bc)
    else
        cg=0.0 
    fi
    index=$(expr $index+1 | bc)
    awk 'NR < 24 { print } ' 0/alpha.water > 0/$i
    sed -i "s/0 0 0 0 0 0 0/0 -3 0 0 1 0 0/g" 0/$i
    sed -i "s/alpha.water/$i/g" 0/$i
    awk 'NR > 23 { print } ' 0/alpha.water > tmp
    awk '/)/ {exit} {print}' tmp > tmp2
    awk '{ print ('$cw'*$l+'$cg'*(1-$l))}' tmp2 >> 0/$i
    rm tmp tmp2
    echo ")" >> 0/$i
    echo ";" >> 0/$i
    tail -27 0_org/Species >> 0/$i
    sed -i "s/cmol/$cw/g" 0/$i
done

cp -r 0_org/U 0/.
sed -i "s/Uinj/$Uinj/g" 0/U

if [ -d "processor0" ]
then
    echo "decomposePar"
    decomposePar -fields > decomposeParRT.out

    rm -rf 0
fi
