#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "process surface area in parallel on $NP processors"
    mpiexec -np $NP processPoroSurf -parallel > processPoroSurf.out
else
    echo -e "process surface area"
    processPoroSurf > processPoroSurf.out
fi
