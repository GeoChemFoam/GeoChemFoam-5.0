#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="4StraightFractures"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'

#Choose if the image is compressed or not
compressed='yes'

# Define image dimensions
x_dim=500
y_dim=500
z_dim=1

#Values of solid, pore, and minimum porosity value for the solid phase (note: this CANNOT be 0)
pores_value=255
solid_value=0
eps_min=0

# tort=1/(eps)^index
tortIndex=1

# Define cropping parameters
x_min=0
x_max=500
y_min=0
y_max=500
z_min=0
z_max=5

#padding for inlet/outlet
padWidth=0

# define resolution (m)
res=0.0001

# number of cells of initial mesh
n_x=250
n_y=250
n_z=1

# flow direction 0 or 1, 2 is empty
direction=0

# Number of processors
NP=8

#### END OF USER INPUT #######################################################

if [ $format != 'raw' ]
then
        echo "ERROR: only raw format is implemented for this solver"
        exit
fi

#Insert dimensions in postProcessDict
x_1=0
y_1=0
z_1=0

xSize=$(echo "$x_max - $x_min" | bc)
ySize=$(echo "$y_max - $y_min" | bc)
zSize=$(echo "$z_max - $z_min" | bc)

x_2=$(expr $xSize*$res | bc)
y_2=$(expr $ySize*$res | bc)
z_2=$(expr $zSize*$res | bc)


cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_1/g" system/postProcessDict
sed -i "s/y_1/$y_1/g" system/postProcessDict
sed -i "s/z_1/$z_1/g" system/postProcessDict

sed -i "s/x_2/$x_2/g" system/postProcessDict
sed -i "s/y_2/$y_2/g" system/postProcessDict
sed -i "s/z_2/$z_2/g" system/postProcessDict

sed -i "s/flowdir/$direction/g" system/postProcessDict

mkdir constant/polyMesh
cp system/blockMeshDict2D$direction system/blockMeshDict

#Dummy fluid properties
cp system/fvSolution1 system/fvSolution

mkdir 0


filename=$Image_name\.$format

if [ $compressed == 'yes' ]
then
        filename=$Image_name\.$format\.tar.gz
fi

mkdir constant/triSurface
cp $dir\/$filename constant/triSurface/.
cd constant/triSurface
if [ $compressed == 'yes' ]
then
        tar -xf $filename
fi
cd ../..


python system/createblockmesh.py --xDim $x_dim --yDim $y_dim --zDim $z_dim --xMin $x_min --xMax $x_max --yMin $y_min --yMax $y_max --zMin $z_min --zMax $z_max --nX $n_x --nY $n_y --nZ $n_z --res $res --Image_name $Image_name --padWidth $padWidth --pores_value $pores_value --solid_value $solid_value --eps_min $eps_min --tortIndex $tortIndex --direction $direction 

rm -rf constant/triSurface

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict
    echo "decomposePar"
    decomposePar > decomposePar.out
    rm -rf 0
fi

echo -e "Case initialised. It is advised to check in paraview to confirm mesh and 0/eps are reasonable before running flow"





