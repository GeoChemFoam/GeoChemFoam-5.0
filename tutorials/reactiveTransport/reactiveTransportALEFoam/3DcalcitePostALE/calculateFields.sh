#!/bin/bash

# Load user environment variables
source $HOME/.bashrc

#export $GCFOAM_DIR/lib

cd ../temp

cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    cp -r 0_org 0

    echo "decompose parallel"
    decomposePar -fields  > decomposePar0.out

    echo "map fields to new mesh in parallel on $NP processors"
    mpiexec -np $NP mapFieldsPar ../3DcalcitePostALE -case ../temp -sourceTime latestTime -parallel > mapFields.out

    for i in processor*; do mv $i/0/pointMotionU* $i/0/pointMotionU; done

    for i in processor*;do cp -r $i/constant/polyMesh $i/0/.; done


    echo "run reactiveTransportALEFoam in parallel on $NP processors for 1e-6 second"
    mpiexec -np $NP reactiveTransportALEFoam -parallel > reactiveTransportALEFoam0.out

    for i in processor*
    do
        cp $i/1e-06/C $i/0/.
        cp $i/1e-06/U $i/0/.
        cp $i/1e-06/p $i/0/.
        cp $i/1e-06/phi $i/0/.
        cp $i/1e-06/pointMotionU $i/0/.
        cp $i/1e-06/cellMotionU $i/0/.
    done
    rm -rf processor*/1e-06
else
    cp -r 0_org 0

    echo "map fields to new mesh"
    mapFieldsGC ../3DcalcitePostALE -case ../temp -sourceTime latestTime > mapFields.out
    mv 0/pointMotionU* 0/pointMotionU

    echo "run reactiveTransportALEFoam for 1e-06 second" 
    reactiveTransportALEFoam > reactiveTransportALEFoam0.out

    cp 1e-06/C 0/.
    cp 1e-06/U 0/.
    cp 1e-06/p 0/.
    cp 1e-06/phi 0/.
    cp 1e-06/pointMotionU 0/.
    cp 1e-06/cellMotionU 0/.

    rm -rf 1e-06

    cp -r constant/polyMesh 0/.

fi
echo -e "Fields have been initialised."
