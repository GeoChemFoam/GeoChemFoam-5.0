#!/bin/bash

# Load user environment variables
source $HOME/.bashrc

#export $GCFOAM_DIR/lib

MPIRUN=mpirun

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
else
    export NP=1
fi

cd ../temp

echo -e "Remove internal faces"
./removeInternalFaces.sh

echo -e "Create stl"
surfaceMeshExtract -latestTime -patches '(movingWalls)' constant/triSurface/Image_meshed.stl > surfaceMeshExtract.out
cp constant/triSurface/Image_meshed.stl ../3DcalcitePostALE/constant/triSurface/Image_meshed.stl
rm -rf 0 0.* *e-* [1-9]*
rm -f *.out
rm -rf processor*
rm -f *.csv
rm -rf constant/polyMesh
rm -rf polyMesh_old

# Create background mesh
echo -e "Create background mesh"
blockMesh  > blockMesh.out

if [ $NP -gt 1 ]
then
    # Decompose background mesh
    echo -e "Decompose background mesh"
    decomposePar > decomposeBlockMesh.out

    # Run snappyHexMesh in parallel
    echo -e "Run snappyHexMesh in parallel on $NP processors"
    $MPIRUN -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out

    # reconstruct mesh to fields decomposition
    echo -e "reconstruct parallel mesh"
    reconstructParMesh -constant > reconstructParMesh.out
else
    echo "Run snappyHexMesh"
    snappyHexMesh -overwrite > snappyHexMesh.out
fi
