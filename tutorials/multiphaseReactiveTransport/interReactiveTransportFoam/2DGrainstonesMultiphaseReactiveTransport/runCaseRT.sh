#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=1.5e-3
WriteTimestep=1e-4
runTimestep=5e-7

#### END OF USER INPUT #######################################################

cp system/fvSolution1 system/fvSolution
cp system/fvSchemes1 system/fvSchemes
cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "Run interReactiveTransportFoam in parallel on $NP processors"
    mpiexec -np $NP interReactiveTransportFoam  -parallel > interReactiveTransportFoamRT.out
else
    echo -e "Run interReactiveTransportFoam"
    interReactiveTransportFoam > interReactiveTransportFoamRT.out
fi

