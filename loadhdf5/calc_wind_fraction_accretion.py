# import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace
import h5py
import ioformat
import sys

# Based on readhdf5.py

# Calculate the fraction of winds from SFRINFO accretion info:
# Use l50n288.
# Halo Gas: Hot/Cold

# 1. Get stars.AccKey from z108
# 2. Find the accretion events (Should do a separate table)
# 3. Put accretion into redshift bins.

NBINS_ASCALE = 100

errormsg = "Usage: calc_wind_fraction_accretion.py modelname lbox"
if(len(sys.argv) != 3):
    raise ValueError(errormsg)
else:
    modelname = sys.argv[1]
    lbox = (float)(sys.argv[2])

def find_abins_from_ascale(a):
    da = 1./ NBINS_ASCALE
    return min((int)(a / da), NBINS_ASCALE-1)

ascales = linspace(0., 1., NBINS_ASCALE+1)
ascales = 0.5 * (ascales[:-1] + ascales[1:])
zred = 1./ascales - 1.

fbase = "/home/shuiyao_umass_edu/scidata/"
# sfrinfobase = fbase + modelname + "/SFRINFO/"
fstarinfo = fbase + modelname + "/" + modelname + "_108.starinfo"

mbins = array([0.1] * NBINS_ASCALE)
mwbins = array([0.0] * NBINS_ASCALE)
mzbins = array([0.0] * NBINS_ASCALE)

tab = genfromtxt(fstarinfo, names=True)

for t in tab:
    bidx = find_abins_from_ascale(t['a_form'])
    if(t['a_last'] > 0): continue
    mbins[bidx] += t['Mass']
    mwbins[bidx] += t['WindMass']
    mzbins[bidx] += t['Mass'] * t['Z']

fout = open("/home/shuiyao_umass_edu/scidata/"+modelname+"/"+modelname+".wfracAcc", "w")
fout.write("#z wind metals\n")
for i in range(NBINS_ASCALE):
    fout.write("%5.3f %7.5e %7.5e\n" % \
               (zred[i], mwbins[i]/mbins[i], mzbins[i]/mbins[i]))
fout.close()    


