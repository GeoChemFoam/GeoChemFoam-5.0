#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "processSpeciesTransfer in parallel on $NP processors"
    mpiexec -np $NP processSpeciesTransfer -parallel > processSpeciesTransferT.out
else
    echo -e "processSpeciesTransfer"
    processSpeciesTransfer > processSpeciesTransferT.out
fi
