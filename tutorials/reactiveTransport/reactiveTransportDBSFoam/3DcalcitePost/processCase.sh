#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "process solid area in parallel on $NP processors"
    mpiexec -np $NP processSolidArea -parallel > processSolidAreaRT.out
else
    echo -e "process solid area"
    processSolidArea > processSolidAreaRT.out
fi
