#!/bin/bash

###### USERS INPUT ############################################################

#Smoothing parameters: smooth surface when image has artifical roughness created ny segmentation to avoid error when using adaptive mesh
nSmooth=1
cSmooth=0.1

#### END OF USER INPUT #######################################################

cp system/fvSolutionInit system/fvSolution
sed -i "s/nSmooth/$nSmooth/g" system/fvSolution
sed -i "s/cSmooth/$cSmooth/g" system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "smoothSolidSurface in parallel on $NP processors"
    mpiexec -np $NP smoothSolidSurface -parallel  > smoothSolidSurface.out
else
    echo -e "smoothSolidSurface"
    smoothSolidSurface > smoothSolidSurface.out
fi




