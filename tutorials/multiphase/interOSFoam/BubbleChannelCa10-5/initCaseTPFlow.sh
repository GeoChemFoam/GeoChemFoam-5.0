#!/bin/bash
###### USERS INPUT ############################################################

#Injection velocity
Ux=0.000167

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    mkdir 0
    cp 0_orig/U 0/.
    sed -i "s/Ux/$Ux/g" 0/U

    echo "decomposePar"
    decomposePar -fields > decomposeParTP.out

    rm -rf 0
else
    cp 0_orig/U 0/.
    sed -i "s/Ux/$Ux/g" 0/U
fi

echo -e "Case initialised"
