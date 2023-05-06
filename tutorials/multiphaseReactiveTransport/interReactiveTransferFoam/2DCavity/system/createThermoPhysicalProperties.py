import numpy as np
import array
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--Species', nargs='+')
parser.add_argument('--Diff', nargs='+')
parser.add_argument('--He', type=float)
parser.add_argument('--Mw', nargs='+')
opt = parser.parse_args()

Species=opt.Species
Diff=opt.Diff
He=opt.He
Mw=opt.Mw

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
f.write("Phase1 phreeqcMixture;"+'\n')
f.write("Phase2 inertMultiComponentMixture;"+"\n")
f.write("PhaseTransfer true;"+"\n")
f.write("solutionSpecies"+'\n')
f.write("{"+'\n')
for i in range (0,nlist):
    f.write('    '+Species[i]+'\n')
    f.write('    {'+'\n')
    f.write('        D1 D1 [0 2 -1 0 0 0 0] '+Diff[i]+';'+'\n')
    f.write('        D2 D2 [0 2 -1 0 0 0 0] 0;'+'\n')
    if (Species[i]=='CO2'):
        f.write('        H H [0 0 0 0 0 0 0] '+str(He)+';'+'\n') 
    else:
        f.write('        H H [0 0 0 0 0 0 0] 1e-9;'+'\n')
    f.write('        Mw Mw [1 0 0 0 -1 0 0] '+Mw[i]+';'+'\n')
    f.write('    }'+'\n')
f.write("}"+'\n')
f.write("surfaceSpecies"+'\n')
f.write("{"+'\n')
f.write("}"+'\n')
f.write("surfaceMasters"+'\n')
f.write("{"+'\n')
f.write("}"+'\n')
f.close()
