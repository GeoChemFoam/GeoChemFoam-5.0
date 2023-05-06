#!/bin/bash

###### USERS INPUT ############################################################

## Define the total time of the simulation and how often to output concentration fields
TotalTime=4
WriteTimestep=0.5
runTimestep=0.005

#### END OF USER INPUT #######################################################

cp system/controlDictReactiveTransport system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # Run reactiveTransportFoam in parallel
    echo -e "Run reactiveTransportFoam in parallel on $NP processors"
    mpirun -np $NP reactiveTransportFoam -parallel  > reactiveTransportFoamRT.out
else
    echo -e "Run reactiveTransportFoam"
    reactiveTransportFoam > reactiveTransportFoamRT.out
fi

