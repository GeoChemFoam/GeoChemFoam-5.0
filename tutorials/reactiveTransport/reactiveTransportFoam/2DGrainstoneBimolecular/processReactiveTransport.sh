#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "Run processConcentration in parallel on $NP processors" 
    mpirun -np $NP processConcentration -parallel > processConcentrationRT.out
else
    echo -e "Run processConcentration" 
    processConcentration > processConcentrationRT.out
fi

