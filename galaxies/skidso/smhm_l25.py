import ioformat
from numpy import histogram, logspace, linspace, array, sqrt, log10
import sys
import os

errormsg = "Usage: "+sys.argv[0]+" modelname"
if(len(sys.argv) != 2):
    print errormsg
    sys.exit(1)
else:
    modelname = sys.argv[1]

DATABASEDIR = "/proj/shuiyao/"
SKIDBASE = modelname+"/"
#SKIDBASE = modelname+"/"+modelname+"/"
#OUTBASE = "/data002/shuiyao/gsmfs_smhms/"+modelname
OUTBASE = "/scratch/shuiyao/sci/PHEW_TEST/"+modelname
BOXSIZE = 25.0
unit_m = 433697.3519675

# z={4.0, 3.0, 1.0, 0.0} -> n={033, 043, 078, 108}
for n in ["033","058","078","098","108"]:
    outfile = OUTBASE+"/smhm_"+n+".txt"

    print outfile

    galfile = DATABASEDIR+SKIDBASE+"gal_z"+n+".stat"
    sofile = DATABASEDIR+SKIDBASE+"so_z"+n+".sovcirc"
    parfile = DATABASEDIR+SKIDBASE+"so_z"+n+".par"

    print "compiled."
    fgal = open(galfile, "r")
    fso = open(sofile, "r")
    fpar = open(parfile, "r")
    f = open(outfile, "w")
    f.write("# Header: Mstar Mvir CentralFlag\n")
    fso.readline()
    for line in fgal:
        Mstar = float(line.split()[4]) * unit_m * 1.e10 / 0.7
        # Mh = 10.**(float(fso.readline().split()[1]))
        spt = fso.readline().split()
        Mh = float(spt[1]) / 0.7
        Msub = float(spt[6]) / 0.7
        # Mh = 10.**float(fso.readline().split()[1]) / 0.7
        spt = fpar.readline().split()
        if(spt[0] == spt[1]):
            f.write(str(Mstar)+" "+str(Mh)+" "+"1 "+str(Msub)+"\n")
        else:
            f.write(str(Mstar)+" "+str(Mh)+" "+"0 "+str(Msub)+"\n")
    f.close()
    fgal.close()
    fso.close()
    fpar.close()
