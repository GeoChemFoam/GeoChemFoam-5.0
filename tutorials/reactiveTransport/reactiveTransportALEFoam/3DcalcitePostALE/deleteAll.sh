#!/bin/bash

set -e

rm -rf 0 0.* *e-* 1* 2* 3* 4* 5* 6* 7* 8* 9*

rm -f *.out
rm -rf processor*
rm -rf constant/polyMesh
rm -f constant/triSurface/Image_meshed*
rm -f constant/triSurface/pore_ind*
rm -rf polyMesh_old
rm -rf 0_org
rm -rf ../temp
rm -f *.csv

