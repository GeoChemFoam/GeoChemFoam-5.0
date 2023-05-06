#!/bin/bash
###### USERS INPUT ############################################################

## postporcessing window
x_1=2e-4
x_2=6e-4
y_1=0
y_2=5e-5

#### END OF USER INPUT #######################################################

cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_1/g" system/postProcessDict
sed -i "s/x_2/$x_2/g" system/postProcessDict
sed -i "s/y_1/$y_1/g" system/postProcessDict
sed -i "s/y_2/$y_2/g" system/postProcessDict


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo "postprocess phase concentration in parallel on $NP processors"
    mpiexec -np $NP processPhaseConcentration -parallel > processPhaseConcentrationRun.out

    echo "postprocess interface transfer in parallel on $NP processors"
    mpiexec -np $NP processInterfaceTransfer -parallel > processInterfaceTransferRun.out
else
    echo "postprocess phase concentration"
    processPhaseConcentration > processPhaseConcentrationRun.out

    echo "postprocess interface transfer"
    processInterfaceTransfer > processInterfaceTransferRun.out
fi


