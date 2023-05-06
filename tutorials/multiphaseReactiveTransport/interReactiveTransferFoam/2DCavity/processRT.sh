#!/bin/bash
###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "process saturation in parallel on $NP processors"
    mpiexec -np $NP processSaturation -parallel > processSatRT.out

    echo "process concentration in parallel on $NP processors"
    mpiexec -np $NP processPhaseConcentration -parallel > processPhaseConcentrationRT.out


else
    echo "process saturation"
    processSaturation > processSatRT.out
    echo "process concentration"
    processPhaseConcentration > processPhaseConcentrationRT.out
fi
