#!/bin/bash

set -e

rm -f *.out
rm -rf constant/polyMesh
rm -rf constant/cell* constant/point*
rm -rf processor*
rm -rf 0.* [1-9]* 
rm -rf 0
rm -rf constant/triSurface/Image_meshed*
rm -f constant/triSurface/pore_ind*
rm -f *.csv
rm -f system/SPDict

