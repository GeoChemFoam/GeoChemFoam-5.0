#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=2000

#### END OF USER INPUT #######################################################

cp system/controlDictFlow system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # Run simpleFoam in parallel
    echo -e "Run simpleFoam in parallel on $NP processors"
    mpirun -np $NP simpleFoam -parallel  > simpleFoamFlow.out
else
    echo -e "Run simpleFoam"
    simpleFoam > simpleFoamFlow.out
fi

echo -e "It is advised to check the simpleFoamFlow.out file to confirm flow has converged before using the permeability output or running transport. If flow has converged the file will say 'SIMPLE solution converged in X iterations'. If flow has not converged, increase TotalTime and/or change p tolerance in sustem/fvSolution and re-run script."

