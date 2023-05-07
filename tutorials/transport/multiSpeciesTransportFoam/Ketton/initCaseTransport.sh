#!/bin/bash

###### USERS INPUT ############################################################

#Number of species
nSpecies=2

##Species name
Species=('A' 'B')


## Define the diffusion coefficient of the species as a solute (m^2/s)
## Eg for water vapor in air this would be 2.42e-5, for CO2(aq) in water this would be 3e-9
Diff=('2.42e-5' '2.42e-7')

#### END OF USER INPUT #######################################################

python system/createThermoPhysicalProperties.py --Species ${Species[@]} --Diff ${Diff[@]}

if [ -d "processor0" ]
then
    mkdir 0
    for i in ${Species[@]}
        do cp 0_orig/Species 0/$i; sed -i "s/Species/$i/" 0/$i
    done
    echo "decomposePar"
    decomposePar -fields > decomposeParTransport.out
    rm -rf 0
else
    for i in ${Species[@]}
        do cp 0_orig/Species 0/$i; sed -i "s/Species/$i/" 0/$i
    done
fi
