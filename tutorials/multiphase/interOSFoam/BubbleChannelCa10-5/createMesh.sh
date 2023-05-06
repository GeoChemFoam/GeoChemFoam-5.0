#!/bin/bash

###### USERS INPUT ############################################################

# Define channel dimensions (m)
x_dim=0.0006
y_dim=0.0001

# number of cells in mesh
n_x=600
n_y=100

#Bubble initially defined as rectangle from x1 to x2
x1=0.00002
x2=0.00022

#number of processores
NP=8

#### END OF USER INPUT #######################################################

MPIRUN=mpirun 

# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict2D system/blockMeshDict

sed -i "s/dx/$x_dim/g" system/blockMeshDict
sed -i "s/dy/$y_dim/g" system/blockMeshDict

sed -i "s/nx/$n_x/g" system/blockMeshDict
sed -i "s/ny/$n_y/g" system/blockMeshDict

cp system/controlDictInit system/controlDict

blockMesh  > blockMesh.out

cp system/setFieldsDict0 system/setFieldsDict1

sed -i "s/ydim/$y_dim/g" system/setFieldsDict1

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict0 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict

    echo "decomposePar"
    decomposePar > decomposePar0.out
fi

echo -e "Mesh Initialised. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
