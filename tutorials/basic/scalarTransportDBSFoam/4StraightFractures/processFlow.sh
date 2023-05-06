#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then 
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    rm -rf processor*/0
    for i in processor*; do mv $i/[1-9]* $i/0; done
    echo -e "processPoroPerm in parallel on $NP processors"
    mpirun -np $NP processPoroPerm -parallel > poroPermFlow.out
else
    rm -rf 0
    mv -f [1-9]* 0
    echo -e "processPoroPerm"
    processPoroPerm > poroPermFlow.out
fi
