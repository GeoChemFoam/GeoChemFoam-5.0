#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "process relative permeability in parallel on $NP processors"
    mpiexec -np $NP processRelPerm -parallel > processRelPermTP.out
else
    echo -e "process relative permeability"
    processRelPerm > processRelPermTP.out
fi


