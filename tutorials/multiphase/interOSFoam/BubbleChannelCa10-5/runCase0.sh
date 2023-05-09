#!/bin/bash
###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=0.00024
runTimestep=4e-5

### OS time-steping parameters ###
#capillary relaxation timestep
dtCR=5e-8
RelaxMin=4
Residual=0.0001

#### END OF USER INPUT #######################################################

cp system/controlDict0 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$TotalTime/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "Run interOSFoam in parallel on $NP processors"
    mpiexec -np $NP interOSFoam  -parallel > interOSFoam0.out

    for i in processor*; do cp $i/$TotalTime/alpha.phase1 $i/0/.; done
    rm -rf processor*/0.*
    rm -rf processor*/[1-9]*
else
    echo "Run interOSFoam"
    interOSFoam > interOSFoam0.out

    cp $TotalTime/alpha.phase1 0/. 

    rm -rf 0.* [1-9]*
fi
