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
    # if PLATFORM is ARCHER2 then use srun, otherwise use mpirun
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
       srun --distribution=block:block --hint=nomultithread -n $NP simpleFoam -parallel  > simpleFoamFlow.out
    else
       mpirun -np $NP simpleFoam -parallel  > simpleFoamFlow.out
    fi
else
    echo -e "Run simpleFoam"
    simpleFoam > simpleFoamFlow.out
fi

echo -e "It is advised to check the simpleFoamFlow.out file to confirm flow has converged before using the permeability output or running transport. If flow has converged the file will say 'SIMPLE solution converged in X iterations'. If flow has not converged, Increase TotalTime and/or change p tolerance in system/fvSolution and re-run the script."

