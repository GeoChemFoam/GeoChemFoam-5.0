import numpy as np
import array
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--Species', nargs='+')
parser.add_argument('--Diff', nargs='+')
opt = parser.parse_args()

Species=opt.Species
Diff=opt.Diff

nlist=len(Species)

print('create thermoPhysicalProperties')
os.system('rm -f constant/thermoPhysicalProperties')
f=open('constant/thermoPhysicalProperties','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       dictionary;"+'\n')
f.write("    object      thermoPhysicalProperties;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("solutionSpecies"+'\n')
f.write("{"+'\n')
for i in range (0,nlist):
    f.write('    '+Species[i]+'\n')
    f.write('    {'+'\n')
    f.write('        D D [0 2 -1 0 0 0 0] '+Diff[i]+';'+'\n')
    f.write('    }'+'\n')
f.write("}"+'\n')
f.close()
