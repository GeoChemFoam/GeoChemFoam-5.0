#!/bin/bash
###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation and how often to output
TotalTime=0.12
initTimestep=0.000001
maxTimestep=0.00004
maxCourant=0.2

#### END OF USER INPUT #######################################################

cp system/controlDict0 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$TotalTime/g" system/controlDict
sed -i "s/initTimestep/$initTimestep/g" system/controlDict
sed -i "s/maxTimestep/$maxTimestep/g" system/controlDict
sed -i "s/maxCourant/$maxCourant/g" system/controlDict

if [ -d "processor0" ]
then 
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "Run interGCFoam in parallel on $NP processors"
    mpiexec -np $NP interGCFoam -parallel > interGCFoam0.out

    for i in processor*; do cp $i/$TotalTime/alpha.phase1 $i/0/.; done
    rm -rf processor*/$TotalTime
else
    echo "Run interGCFoam"
    interGCFoam > interGCFoam0.out

    cp $TotalTime/alpha.phase1 0/. 
    rm -rf $TotalTime
fi

