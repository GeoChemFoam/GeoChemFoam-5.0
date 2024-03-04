#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "Run snappyHexMesh in parallel on $NP processors"
# if PLATFORM is ARCHER2 then use srun, otherwise use mpirun
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
       srun --distribution=block:block --hint=nomultithread -n $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out
    else
       mpirun -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out
    fi
else
    echo "Run snappyHexMesh"
    snappyHexMesh -overwrite > snappyHexMesh.out
fi

echo -e "Mesh created. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
