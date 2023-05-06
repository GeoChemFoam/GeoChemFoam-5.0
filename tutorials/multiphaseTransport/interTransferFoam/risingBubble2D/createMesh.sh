#!/bin/bash
###### USERS INPUT ############################################################

# Define channel dimensions (m)
x_dim=0.012
y_dim=0.024

# number of cells in mesh
n_x=120
n_y=240

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

cp system/controlDictInit system/controlDict

blockMesh  > blockMesh.out

if [ $NP -gt 1 ]
then
    echo "decomposePar"
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict
    decomposePar > decomposePar.out
fi

echo -e "Mesh created. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
