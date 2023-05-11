#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "Run snappyHexMesh in parallel on $NP processors"
    mpirun -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out
else
    echo "Run snappyHexMesh"
    snappyHexMesh -overwrite > snappyHexMesh.out
fi

echo -e "Mesh created. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
