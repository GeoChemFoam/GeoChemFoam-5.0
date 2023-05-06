#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=50000
WriteTimestep=500
runTimestep=5

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

cp system/fvSolution2 system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "Run scalarTransportDBSFoam in parallel on $NP processor"
    mpiexec -np $NP scalarTransportDBSFoam -parallel  > scalarTransportDBSFoamT.out 
else
    echo -e "Run scalarTransportDBSFoam"
    scalarTransportDBSFoam > scalarTransportDBSFoamT.out
fi
