#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="HM12_2_4"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'
#format='stl'

compressed='yes'

#If format='raw', Define image dimensions
x_dim=2000
y_dim=2000
z_dim=10

#If format='raw', define value of the pores and solid in the image
pores_value=255
solid_value=0

#if format= 'stl' enter the location of a pore
#pore_index_X=0
#pore_index_Y=0
#pore_index_Z=0

# define resolution (m)
res=0.00003

# Define cropping parameters
x_min=700
x_max=1300
y_min=700
y_max=1300
z_min=1
z_max=6

# number of cells of initial mesh
# n*(nlevel+1) should be equal to image dimension when not binning 
n_x=300
n_y=300
#In 2D images, z has to be the empty direction
n_z=1

# flow direction 0 or 1, 2 is empty
direction=0

# Percent of image cut for permeability calculation (boundary and capillary end effect)
cut_x=0.25
cut_y=0.05

# Number of processors
NP=8

#### END OF USER INPUT #######################################################

MPIRUN=mpirun 

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
    surfaceTransformPoints -translate '(-0.5 -0.5 -0.5)' Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints.out
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
    pore_index_2=$(expr $res*$pore_index_Z-$res*$x_min | bc)
fi

echo -e "Coordinates at center of a pore = ($pore_index_0,$pore_index_1,$pore_index_2)" 

cd ../..
 
# Create background mesh
echo -e "Create background mesh"
cp system/blockMeshDict2D$direction system/blockMeshDict
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

cp system/controlDictInit system/controlDict

blockMesh  > blockMesh.out

x_1=$(expr $cut_x*$dx*$res | bc)
y_1=$(expr $cut_y*$dy*$res | bc)
z_1=0

x_2=$(expr $dx*$res-$x_1 | bc)
y_2=$(expr $dy*$res-$y_1 | bc)
z_2=$(expr $dz*$res | bc)


cp system/postProcessDict1 system/postProcessDictRun
sed -i "s/x_1/$x_1/g" system/postProcessDictRun
sed -i "s/y_1/$y_1/g" system/postProcessDictRun
sed -i "s/z_1/$z_1/g" system/postProcessDictRun

sed -i "s/x_2/$x_2/g" system/postProcessDictRun
sed -i "s/y_2/$y_2/g" system/postProcessDictRun
sed -i "s/z_2/$z_2/g" system/postProcessDictRun

sed -i "s/flowdir/$direction/g" system/postProcessDictRun


cp system/decomposeParDictRun system/decomposeParDict
sed -i "s/NP/$NP/g" system/decomposeParDict


# Run snappyHexMesh in parallel
cp system/snappyHexMeshDict1 system/snappyHexMeshDict

sed -i "s/poreIndex0/$pore_index_0/g" system/snappyHexMeshDict
sed -i "s/poreIndex1/$pore_index_1/g" system/snappyHexMeshDict
sed -i "s/poreIndex2/$pore_index_2/g" system/snappyHexMeshDict

xmax=$(expr $dx*$res | bc)
ymax=$(expr $dy*$res | bc)
zmax=$(expr $dz*$res | bc)

if [ $direction -eq 0 ]
then
    linter=$x_1
elif [ $direction -eq 1 ]
then
    linter=$y_1
fi

echo "$x_max $y_max $z_max $n_x $n_y $n_z $direction $linter" > system/dimensions

if [ $NP -gt 1 ]
then
    # Decompose background mesh
    echo -e "Decompose background mesh"
    decomposePar > decomposeBlockMesh.out

    rm -rf constant/polyMesh
fi

echo -e "BlockMesh and image created. Ready for snappyHexMesh" 
