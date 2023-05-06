#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "processPhaseConcentration in parallel on $NP processors"
    mpiexec -np $NP processPhaseConcentration -parallel > processPhaseConcentrationRT.out
else
    echo -e "processPhaseConcentration"
    processPhaseConcentration > processPhaseConcentrationRT.out
fi

