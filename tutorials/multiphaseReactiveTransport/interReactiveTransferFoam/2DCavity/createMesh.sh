#!/bin/bash

###### USERS INPUT ############################################################

# Define channel dimensions (m)
x_dim=0.006
y_dim=0.001
z_dim=0.001

# number of cells in channel
n_x=60
n_y=10

## Define cavity dimensions
## The cavity will be meshed with the same resolution as the channel
xc_dim=0.002
yc_dim=0.002

## Number of processors
NP=4

#### END OF USER INPUT #######################################################

lx1=$(expr 0.5*$x_dim-0.5*$xc_dim | bc -l)
nx1=$(expr $lx1*$n_x/$x_dim | bc)
nx1=${nx1%.*}

lx2=$(expr $lx1+$xc_dim | bc -l)
nx2=$(expr $xc_dim*$n_x/$x_dim | bc)
nx2=${nx2%.*}

lx3=$x_dim
nx3=$(expr $n_x-$nx1-$nx2 | bc)

ly1=$y_dim
ny1=$n_y

ly2=$(expr -1*$yc_dim | bc -l)
ny2=$(expr $yc_dim/$y_dim*$n_y | bc)
ny2=${ny2%.*}

ly3=$(expr -1*$x_dim/$n_x | bc -l)


# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict2D system/blockMeshDict

sed -i "s/lx1/$lx1/g" system/blockMeshDict
sed -i "s/lx2/$lx2/g" system/blockMeshDict
sed -i "s/lx3/$lx3/g" system/blockMeshDict
sed -i "s/ly1/$ly1/g" system/blockMeshDict
sed -i "s/ly2/$ly2/g" system/blockMeshDict
sed -i "s/lz/$z_dim/g" system/blockMeshDict


sed -i "s/nx1/$nx1/g" system/blockMeshDict
sed -i "s/nx2/$nx2/g" system/blockMeshDict
sed -i "s/nx3/$nx3/g" system/blockMeshDict
sed -i "s/ny1/$ny1/g" system/blockMeshDict
sed -i "s/ny2/$ny2/g" system/blockMeshDict


cp system/controlDictInit system/controlDict

echo "blockMesh"
blockMesh  > blockMesh.out

echo "$lx1 $lx2 $ly2 $ly3 0 $z_dim" > system/dimensions

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict

    echo "decomposePar"
    decomposePar > decomposePar.out
fi

echo -e "Mesh created. It is advised to check in paraview to confirm mesh of porespace is reasonable before running flow" 
