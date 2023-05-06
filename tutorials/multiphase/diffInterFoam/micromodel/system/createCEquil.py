import numpy as np
import array
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--xMin', type=float, help='xMin')
parser.add_argument('--xMax', type=float, help='xMax')
parser.add_argument('--yMin', type=float, help='yMin')
parser.add_argument('--yMax', type=float, help='yMax')
parser.add_argument('--zMin', type=float, help='zMin')
parser.add_argument('--zMax', type=float, help='zMax')
parser.add_argument('--direction', type=float, help='direction')
parser.add_argument('--linter',type=float, help='linter')
parser.add_argument('--nX', type=int, help='nX')
parser.add_argument('--nY', type=int, help='nY')
parser.add_argument('--nZ', type=int, help='nZ')
parser.add_argument('--eps', type=float, help='eps')
opt = parser.parse_args()

xMin=opt.xMin
xMax=opt.xMax
yMin=opt.yMin
yMax=opt.yMax
zMin=opt.zMin
zMax=opt.zMax

direction=opt.direction
linter=opt.linter

nX=opt.nX
nY=opt.nY
nZ=opt.nZ

eps=opt.eps

os.system('echo "process mesh center"')
os.system('processMeshCellCenters > processMeshCellCenter.out')

x = np.zeros(nX*nY*nZ,dtype=float)
y = np.zeros(nX*nY*nZ,dtype=float)
z = np.zeros(nX*nY*nZ,dtype=float)


file = open("0/cellCenters","r")
Lines = file.readlines()
count =0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      x[count]=float(ls.split("(")[1].split(")")[0].split()[0])
      y[count]=float(ls.split("(")[1].split(")")[0].split()[1])
      z[count]=float(ls.split("(")[1].split(")")[0].split()[2])
      count +=1
  if (ls=="("):
      wbool=1

ncell = count

os.system('rm -f 0/cellCenters')


C=np.zeros(ncell,dtype=float)
print('calculate C')
for i in range (0,ncell):
    if (direction==0):
        C[i] = np.tanh((x[i]-linter)/np.sqrt(2)/eps) 
    if (direction==1):
        C[i] = np.tanh((y[i]-linter)/np.sqrt(2)/eps)
    if (direction==2):
        C[i] = np.tanh((z[i]-linter)/np.sqrt(2)/eps)


print('create C')
f=open('0/C','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       volScalarField;"+'\n')
f.write("    object      C;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("dimensions      [0 0 0 0 0 0 0];"+'\n')
f.write("internalField   nonuniform List<scalar>"+'\n')
f.write(str(ncell)+'\n')
f.write("("+'\n')
for i in range (0,ncell):
        f.write(str(C[i])+'\n')
f.write(")"+'\n')
f.write(";"+'\n')
f.write(""+'\n')
f.write("boundaryField"+'\n')
f.write("{"+'\n')
f.write("    walls"+'\n')
f.write("    {"+'\n')
f.write("        type wallEnergyConstantContactAngle;"+'\n')
f.write("        theta0 contactAngle;"+'\n')
f.write("    }"+'\n')
f.write("    solidwalls"+'\n')
f.write("    {"+'\n')
f.write("        type wallEnergyConstantContactAngle;"+'\n')
f.write("        theta0 contactAngle;"+'\n')
f.write("    }"+'\n')
f.write("    inlet"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    outlet"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    frontAndBack"+'\n')
f.write("    {"+'\n')
f.write("        type empty;"+'\n')
f.write("    }"+'\n')
f.write("}"+'\n')
f.close()
