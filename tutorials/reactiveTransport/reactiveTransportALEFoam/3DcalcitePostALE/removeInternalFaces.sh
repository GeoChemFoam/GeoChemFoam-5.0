#!/bin/bash




# ###### DO NOT MAKE CHANGES FROM HERE ###################################

set -e

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    mpiexec -np $NP topoSet -parallel > toposet.out
    mpiexec -np $NP createPatch -overwrite -parallel > createPatch.out
    reconstructParMesh -constant > reconstructPar.out
else
    topoSet > toposet.out
    createPatch -overwrite > createPatch.out
fi    


