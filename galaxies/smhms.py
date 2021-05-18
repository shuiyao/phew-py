# OPTIMIZED FOR SMF

import ioformat
import matplotlib.pyplot as plt
from numpy import histogram, logspace, linspace, array, sqrt, log10
from myinit import *

# Compare SMF at different redshift.
# Write smf_###.txt file

# modelname = "l50n576-phew-m5"
modelname = "l50n288-phew-m5-sph"
DATABASEDIR = DIRS['DATA']
SKIDBASE = DIRS['DATA'] + modelname + "/"
# BOXSIZE = 18.0
# unit_m = 3469578.81574 * (18./50.)**3
BOXSIZE = 50.0
unit_m = 3469578.81574

def batch():
    for n in ["033", "058", "078", "108", "098"]:
        galfile = SKIDBASE+"gal_z"+n+".stat"
        sofile = SKIDBASE+"so_z"+n+".sovcirc"
        parfile = SKIDBASE+"so_z"+n+".par"
        outfile = DIRS['SCIDATA']+modelname+"/smhm_"+n+".txt"    

        print ("read: ", sofile)
        print ("write: ", outfile)

        fgal = open(galfile, "r")
        fso = open(sofile, "r")
        fpar = open(parfile, "r")
        f = open(outfile, "w")
        f.write("# Header: Mstar Mvir CentralFlag\n")
        fso.readline()
        for line in fgal:
            Mstar = float(line.split()[4]) * unit_m * 1.e10 / 0.7
            Mh = float(fso.readline().split()[1]) / 0.7
            # Mh = 10.**float(fso.readline().split()[1]) / 0.7
            spt = fpar.readline().split()
            if(spt[0] == spt[1]):
                f.write(str(Mstar)+" "+str(Mh)+" "+"1"+"\n")
            else:
                f.write(str(Mstar)+" "+str(Mh)+" "+"0"+"\n")
        f.close()
        fgal.close()
        fso.close()
        fpar.close()
