#!/bin/bash

###### USERS INPUT ############################################################

TotalTime=800
WriteTimestep=800
initTimestep=1
maxTimestep=50

#### END OF USER INPUT #######################################################

cp system/controlDictRun system/controlDict
cp system/fvSolutionRun system/fvSolution

sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/initTimestep/$initTimestep/g" system/controlDict
sed -i "s/maxTimestep/$maxTimestep/g" system/controlDict


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "run reactiveTransportDBSFoam in parallel on $NP processors"
    mpiexec -np $NP reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoamRT.out
else
    echo -e "run reactiveTransportDBSFoam"
    reactiveTransportDBSFoam > reactiveTransportDBSFoamRT.out 
fi

