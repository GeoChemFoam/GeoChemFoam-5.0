#!/bin/bash

###### USERS INPUT ############################################################

## Define the total time of the simulation and how often to output fields
## Define initial and maximum time-step
TotalTime=800
WriteTimestep=800
initTimestep=1
maxTimestep=50

#### END OF USER INPUT #######################################################

CWD=$(pwd)
set -e

cp system/fvSolutionRun system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
else
    export NP=1
fi

rm -rf polyMesh_old
rm -rf processor*/polyMesh_old

cp -r constant/polyMesh polyMesh_old

if [ -d "processor0" ]
then
    for i in processor*;do cp -r $i/constant/polyMesh $i/polyMesh_old; done
fi

rm -rf ../temp
cp -r $CWD ../temp

cd $CWD

python << END
import os

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


a=0
s=str(a)
removeOld = 0

dirName=os.getcwd()
if $NP>1:
    dirName=dirName+'/processor0'

for directories in os.listdir(dirName):
  if (is_number(directories)):
    if (float(directories)>a):
      a=float(directories)
      s=directories

b=a

if $NP>1:
    for i in range(0,$NP):
        os.system('cp processor'+str(i)+'/'+s+'/polyMesh/* processor'+str(i)+'/polyMesh_old/.')
else:
    os.system('cp '+s+'/polyMesh/* polyMesh_old/.')

print('Time='+str(s)+' seconds ###############################################################################')
while b<round($TotalTime,13):
  b=round(b+$WriteTimestep,13)
  print('Next report step at T='+str(b)+' seconds')
  os.system('cp system/controlDictRun system/controlDict')
  os.system('sed -i "s/var/'+str(b)+'/g" system/controlDict') 
  os.system('sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict') 
  os.system('sed -i "s/initTimestep/$initTimestep/g" system/controlDict') 
  os.system('sed -i "s/maxTimestep/$maxTimestep/g" system/controlDict') 
  while a<b:
    if $NP>1:
        print("run reactiveTransportALEFoam in parallel on $NP processors")
        os.system('mpiexec -np $NP reactiveTransportALEFoam -parallel > reactiveTransportALEFoamRT.out')
    else:
        print("run reactiveTransportALEFoam")
        os.system('reactiveTransportALEFoam > reactiveTransportALEFoamRT.out')
    if removeOld == 1:
        os.system('rm -rf processor*/'+s)
        os.system('rm -rf '+s)
    dirName=os.getcwd()
    if $NP>1:
        dirName=dirName+'/processor0'
    for directories in os.listdir(dirName): 
      if (is_number(directories)):
        if (float(directories)>a):
          a=float(directories)
          s=directories
    if $NP>1:
        for i in range(0,$NP):
            os.system('cp -r processor'+str(i)+'/'+s+'/polyMesh/points processor'+str(i)+'/polyMesh_old/.')
            os.system('cp -r processor'+str(i)+'/polyMesh_old/* processor'+str(i)+'/'+s+'/polyMesh/.')
            os.system('cp -r processor'+str(i)+'/'+s+' ../temp/processor'+str(i)+'/.')
    else:
        os.system('cp '+s+'/polyMesh/points polyMesh_old/.')
        os.system('cp polyMesh_old/* '+s+'/polyMesh/.')
        os.system('cp -r '+s+' ../temp/.')
    print('Remesh')
    os.system( './remesh.sh') 
    os.system('./calculateFields.sh')
    os.system('cd ../3DcalcitePostALE')
    os.system('rm -rf processor*/'+s)
    os.system('rm -rf '+s)
    if $NP>1:
        for i in range(0,$NP):
            os.system('mv ../temp/processor'+str(i)+'/0 processor'+str(i)+'/'+s) 
            os.system('cp processor'+str(i)+'/'+s+'/polyMesh/* processor'+str(i)+'/polyMesh_old/.')
    else:
        os.system('mv ../temp/0 ./'+s)
        os.system('cp '+s+'/polyMesh/* polyMesh_old/.')
    print(' ')
    print('Time='+str(s)+' seconds ###############################################################################')
    print('Next report step at T='+str(b)+' seconds')
    if a<b:
      removeOld=1
    else:
      removeOld=0
  print('Report time '+str(b)+' seconds completed')
END
