#!/bin/bash

###### USERS INPUT ############################################################

#Name of your image (image must be in raw format 8 bit)
Image_name="Est_3phase500cubed4micron"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'

#Choose if the image is compressed or not
compressed='yes'

# Define image dimensions
x_dim=500
y_dim=500
z_dim=500

#Values of solid, pore, and minimum porosity value for the solid phase (note: this CANNOT be 0) 
pores_value=1
solid_value=3

#define the labels of the phases
phases=(1 2 3)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=('1' '0.35' '0.0001')

#define the permeability of each label (note: solid phase should be < 1e-20, pore should be > 1e6)
micro_k=('1e13' '1.82e-15' '1e-20')

# Define cropping parameters
x_min=100
x_max=300
y_min=100
y_max=300
z_min=100
z_max=300

#padding for inlet/outlet (minimum 2)
padWidth=4

# define resolution (m)
res=0.000004

# number of cells of initial mesh
n_x=100
n_y=100
n_z=100

#Mesh refinement level
nlevel=1
#0=no refinment in pores, 1=refinemnet in pores
refineStokes=1

#Choose the direction of flow, 0 is x, 1 is y, 2 is z
direction=0

##Number of processors
NP=4

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
cp system/blockMeshDict$direction system/blockMeshDict

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/1/g" system/controlDict

cp constant/fvOptions$direction constant/fvOptionsRun

cp system/fvSolution1 system/fvSolution
sed -i "s/nSmooth/1/g" system/fvSolution
sed -i "s/cSmooth/0.5/g" system/fvSolution

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

mkdir 0

python system/createblockmesh.py --xDim $x_dim --yDim $y_dim --zDim $z_dim --xMin $x_min --xMax $x_max --yMin $y_min --yMax $y_max --zMin $z_min --zMax $z_max --nX $n_x --nY $n_y --nZ $n_z --nLevel $nlevel --refineStokes $refineStokes --res $res --Image_name $Image_name --padWidth $padWidth --pores_value $pores_value --solid_value $solid_value --direction $direction --micro_por ${micro_por[@]} --micro_k ${micro_k[@]} --phases ${phases[@]}

if [ $NP -gt 1 ]
then
    cp system/decomposeParDict1 system/decomposeParDict
    sed -i "s/NP/$NP/g" system/decomposeParDict
    echo -e "decomposePar"
    decomposePar > decomposePar.out

    rm -rf 0
fi


rm -rf constant/triSurface

echo -e "Mesh created. It is advised to check in paraview to confirm mesh and 0/eps are reasonable before running flow"

