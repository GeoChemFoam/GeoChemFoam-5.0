#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=2000

#### END OF USER INPUT

cp system/fvSolutionSP system/fvSolution
cp system/fvSchemesSP system/fvSchemes
cp system/controlDictSP system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict


if [ -d "processor0"  ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # Run simpleFoam in parallel
    echo -e "Run simpleFoam in parallel on $NP processors"
    mpiexec -np $NP simpleFoam -parallel  > simpleFoamSP.out
else
	echo -e "Run simpleFoam"
	simpleFoam > simpleFoamSP.out
fi

echo -e "Note: Please check the last line of simpleFoamSP.out to confirm the flow field has converged. If it has not, change TotalTime and/or the p tolerance and residual controls in system/fvSolutionSP and rerun the script"



