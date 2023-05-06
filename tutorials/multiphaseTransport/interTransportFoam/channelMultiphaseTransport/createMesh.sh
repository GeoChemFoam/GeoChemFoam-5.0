#!/bin/bash

###### USERS INPUT ############################################################

# Define channel dimensions (m)
x_dim=0.0008
y_dim=0.00005

# number of cells in mesh
n_x=200
n_y=16

#Grading='(1 1 1)'
Grading='(32 50 4)(68 50 1)'
#Grading='(16 25 4)(68 50 1)(16 25 0.25)'

#number of processores
NP=4

#### END OF USER INPUT #######################################################

# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict2D system/blockMeshDict

sed -i "s/dx/$x_dim/g" system/blockMeshDict
sed -i "s/dy/$y_dim/g" system/blockMeshDict

sed -i "s/nx/$n_x/g" system/blockMeshDict
sed -i "s/ny/$n_y/g" system/blockMeshDict

sed -i "s/grading/$Grading/g" system/blockMeshDict

cp system/controlDictInit system/controlDict

blockMesh  > blockMesh.out

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict

    echo -e "decomposePar"
    decomposePar > decomposePar.out
fi

echo -e "Mesh created. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
