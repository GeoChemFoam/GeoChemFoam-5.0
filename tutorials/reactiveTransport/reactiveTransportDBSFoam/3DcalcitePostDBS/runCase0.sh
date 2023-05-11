#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "run reactiveTransportDBSFoam for 1e-06 sec in parallel on $NP processors"
    mpiexec -np $NP reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoam0.out

    echo -e "move 1e-06 to 0"
    for i in processor*; do cp "$i/1e-06/U" "$i/0/."; done
    for i in processor*; do cp "$i/1e-06/p" "$i/0/."; done
    for i in processor*; do cp "$i/1e-06/phi" "$i/0/."; done
    for i in processor*; do cp "$i/1e-06/C" "$i/0/."; done
    for i in processor*; do cp "$i/1e-06/R" "$i/0/."; done
    rm -rf processor*/1e-06
else
    echo -e "run reactiveTransportDBSFoam for 1e-06"
    reactiveTransportDBSFoam > reactiveTransportDBSFoam0.out

    cp -f 1e-06/C 0/.
    cp -f 1e-06/U 0/.
    cp -f 1e-06/p 0/.
    cp -f 1e-06/phi 0/.
    cp -f 1e-06/R 0/.
    rm -rf 1e-06
fi

echo -e "Note: Please check the last line of reactiveTransportDBSFoam0.out to confirm the equation have converged. If it has not, re-run script and/or change the tolerance and residual controls in system/fvSolution" 


