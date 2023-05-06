#!/bin/bash
###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "processSaturation in parallel on $NP processors"
    mpiexec -np $NP processSaturation -parallel > processSatT.out
else
    echo "processSaturation"
    processSaturation > processSatT.out
fi




