#!/bin/bash
###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

grep -oP 'finished in \s*\K\d+' interOSFoamTP.out > relaxSteps.csv

