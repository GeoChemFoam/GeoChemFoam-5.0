#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="calcitePost"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'

#Choose if the image is compressed or not
compressed='yes'

#If format='raw', define image dimensions
x_dim=536
y_dim=300
z_dim=40

#If format='raw', define value of the pores and solid in the image
pores_value=255
solid_value=0

#if format= 'stl' enter the location of a pore
#pore_index_X=0.0001522 #0
#pore_index_Y=0.0001522 #0
#pore_index_Z=0

# define resolution (m)
res=0.000005

# Define cropping parameters
x_min=0
x_max=536
y_min=0
y_max=300
z_min=0
z_max=40

# number of cells of initial mesh
n_x=134
n_y=75
n_z=10

# Level of refinement - mesh is refined by n_level at the pore surfaces
n_level=1

# flow direction (0 is x, 1 is y and 2 is z)
direction=0

# Number of processors
NP=8

#### END OF USER INPUT #######################################################

if [ $format != 'raw' ] && [ $format != 'stl' ]
then
        echo "ERROR: only raw and stl format are implemented for this solver"
        exit
fi

filename=$Image_name\.$format

if [ $compressed == 'yes' ]
then
        filename=$Image_name\.$format\.tar.gz
fi

cp $dir\/$filename constant/triSurface/.
cd constant/triSurface
if [ $compressed == 'yes' ]
then
        tar -xf $filename
fi

if [ $format == 'raw' ]
then
    echo -e "make stl"
    python raw2stl.py --x_min=$x_min --x_max=$x_max --y_min=$y_min --y_max=$y_max --z_min=$z_min --z_max=$z_max --pores_value=$pores_value --solid_value=$solid_value  --image_name=$Image_name --x_dim=$x_dim --y_dim=$y_dim --z_dim=$z_dim
    rm $Image_name\.*

    surfaceTransformPoints -translate '(-0.5 -0.5 -0.5)' Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints1.out
    vector="($res $res $res)"
    surfaceTransformPoints -scale "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints2.out
    pore_index_X="$(cat pore_indx)"
    pore_index_Y="$(cat pore_indy)"
    pore_index_Z="$(cat pore_indz)"

    pore_index_0=$(expr $pore_index_X*$res | bc)
    pore_index_1=$(expr $pore_index_Y*$res | bc)
    pore_index_2=$(expr $pore_index_Z*$res | bc)


elif [ $format == 'stl' ]
then
    echo -e "reading stl"
    mv $Image_name\.stl Image_meshed.stl
    rm -f $Image_name\.*
    tranX=$(expr -1*$x_min | bc)
    tranY=$(expr -1*$y_min | bc)
    tranZ=$(expr -1*$z_min | bc)
    vector="($tranX $tranY $tranZ)"
    surfaceTransformPoints -translate "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints1.out
    vector="($res $res $res)"
    surfaceTransformPoints -scale "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints2.out

    pore_index_0=$(expr $res*$pore_index_X-$res*$x_min | bc)
    pore_index_1=$(expr $res*$pore_index_Y-$res*$y_min | bc)
    pore_index_2=$(expr $res*$pore_index_Z-$res*$z_min | bc)
fi

echo -e "Coordinates at center of a pore = ($pore_index_0,$pore_index_1,$pore_index_2) in the cropped image"

cd ../..

 
# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict$direction system/blockMeshDict

dx=$(expr $x_max-1*$x_min | bc)
dy=$(expr $y_max-1*$y_min | bc)
dz=$(expr $z_max-1*$z_min | bc)

sed -i "s/dx/$dx/g" system/blockMeshDict
sed -i "s/dy/$dy/g" system/blockMeshDict
sed -i "s/dz/$dz/g" system/blockMeshDict

sed -i "s/nx/$n_x/g" system/blockMeshDict
sed -i "s/ny/$n_y/g" system/blockMeshDict
sed -i "s/nz/$n_z/g" system/blockMeshDict

sed -i "s/res/$res/g" system/blockMeshDict

cp system/controlDict0 system/controlDict

blockMesh  > blockMesh.out

cp system/snappyHexMeshDict1 system/snappyHexMeshDict

sed -i "s/nlevel/$n_level/g" system/snappyHexMeshDict

sed -i "s/poreIndex0/$pore_index_0/g" system/snappyHexMeshDict
sed -i "s/poreIndex1/$pore_index_1/g" system/snappyHexMeshDict
sed -i "s/poreIndex2/$pore_index_2/g" system/snappyHexMeshDict

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

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict

    # Decompose background mesh
    echo -e "Decompose background mesh"
    decomposePar > decomposeBlockMesh.out

    rm -rf constant/polyMesh
fi

echo -e "BlockMesh and Image created. Ready for snappyHeMesh" 
