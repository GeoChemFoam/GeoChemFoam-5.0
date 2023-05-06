#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=2000

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    # Run simpleFoam in parallel
    echo -e "Run laplacianGCFoam in parallel on $NP processors"
    mpirun -np $NP laplacianGCFoam -parallel  > laplacianGCFoamT.out
else
    echo -e "Run laplacianGCFoam"
    laplacianGCFoam > laplacianGCFoamT.out
fi

echo -e "It is advised to check the laplacianGCFoamT.out file to confirm flow has converged. If it has converged the file will say 'SIMPLE solution converged in X iterations'. If it has not converged, Increase TotalTime and/or change T tolerance in system/fvSolution and re-run the script."

