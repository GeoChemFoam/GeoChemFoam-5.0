#!/bin/bash
###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=0.01
WriteTimestep=0.01
runTimestep=0.00004

### OS time-steping parameters ###
#capillary relaxation timestep
dtCR=5e-8
RelaxMin=15
Residual=0.0001

#### END OF USER INPUT #######################################################

cp system/controlDict0 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

cp system/fvSolution0 system/fvSolution
sed -i "s/dtCR/$dtCR/g" system/fvSolution
sed -i "s/RelaxMin/$RelaxMin/g" system/fvSolution
sed -i "s/Residual/$Residual/g" system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "Run interOSFoam in parallel on $NP processors"
    mpiexec -np $NP interOSFoam  -parallel > interOSFoamTP.out
else
   echo "Run interOSFoam"
   interOSFoam > interOSFoamTP.out
fi
