#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

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
    echo "run reactiveTransportALEFoam for 1e-6 second"
    reactiveTransportALEFoam > reactiveTransportALEFoam0.out

    cp 1e-06/C 0/.
    cp 1e-06/U 0/.
    cp 1e-06/p 0/.
    cp 1e-06/phi 0/.
    cp 1e-06/pointMotionU 0/.
    cp 1e-06/cellMotionU 0/.

    rm -rf 1e-06 
fi

echo -e "Note: Please check the last line of reactiveTransportALEFoam0.out to confirm the equation have converged. If it has not, re-run script and/or change the tolerance and residual controls in system/fvSolution" 
