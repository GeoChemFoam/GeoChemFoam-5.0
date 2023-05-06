#!/bin/bash
###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=0.2
WriteTimestep=0.01
initTimestep=5e-7
maxTimestep=1e-5

#### END OF USER INPUT #######################################################

set -e

cp system/controlDict1 system/controlDict

sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/initTimestep/$initTimestep/g" system/controlDict
sed -i "s/maxTimestep/$maxTimestep/g" system/controlDict


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "Run interReactiveTransferFoam in parallel on $NP processors"
    mpiexec -np $NP interReactiveTransferFoam -parallel > interReactiveTransferFoamRT.out
else
    echo "Run interReactiveTransferFoam in parallel on $NP processors"
    interReactiveTransferFoam > interReactiveTransferFoamRT.out
fi
