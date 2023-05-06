import numpy as np
import array
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--Species', nargs='+')
parser.add_argument('--Cmol', nargs='+')
parser.add_argument('--Mw', nargs='+')
parser.add_argument('--rho2',type=float)
opt = parser.parse_args()

Species=opt.Species
Cmol=opt.Cmol
Mw=opt.Mw
rho2=opt.rho2

nlist=len(Species)

print('create setFieldsDict')
os.system('rm -f system/setFieldsDict')
f=open('system/setFieldsDict','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       dictionary;"+'\n')
f.write("    object      setFieldsDict;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("defaultFieldValues"+'\n')
f.write("("+'\n')
f.write('    volScalarFieldValue alpha.water 1'+'\n')
for i in range (0,nlist):
    f.write('    volScalarFieldValue '+Species[i]+' '+Cmol[i]+'\n')
f.write(");"+'\n')
f.write("regions"+'\n')
f.write("("+'\n')
f.write("    boxToCell"+'\n')
f.write("    {"+'\n')
f.write("        box (lx1 ly1 lz1) (lx2 ly2 lz2);"+'\n')
f.write("        fieldValues"+'\n')
f.write("        ("+'\n')
f.write("            volScalarFieldValue alpha.water 0"+'\n')
for i in range (0,nlist):
    c=0
    if (Species[i]=='CO2'):
        c=rho2/float(Mw[i])
    f.write("            volScalarFieldValue "+Species[i]+" "+str(c)+'\n')
f.write("        );"+'\n')
f.write("    }"+'\n')
f.write(");")
f.close()
