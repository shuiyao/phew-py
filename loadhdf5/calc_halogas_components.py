# import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace
import h5py
import ioformat
import sys
# from scipy.interpolate import interp1d
from numpy import polyfit, polyval

modelname = 'l50n288-phewoff'
# Based on readhdf5.py
# if(len(sys.argv) <= 1):
#     raise ValueError("Need to specify modelname.")
# else:
#     modelname = sys.argv[1]

if(len(sys.argv) == 3):
    snapstr = ("000"+sys.argv[2])[-3:]
    print("Set snapstr = ", snapstr)
else:
    print("Use default snapstr = 108")
    snapstr = "108"

BOXSIZE = 50000
UNIT_MASS = 3469581.88
HPARAM = 0.7
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
UNIT_T = 1.e10 * (2./3.) * 1.6726e-24 / 1.3806e-16
TMAX = 5.0
suffix = ""
# TMAX = 4.9
# suffix = ".T80k"

# UNIT_MASS = 3469581.88 * (12./50.) ** 3

#modelname = "l50n288-phewoff"
#modelname = "l50n288-phew-m5"
#modelname = "p50n288fiducial"
fbase = "/nas/astro-th/shuiyao/"
snapname = fbase+modelname+"/snapshot_"+snapstr+".hdf5"
galname = fbase+modelname+"/gal_z"+snapstr+".stat"
soname = fbase+modelname+"/so_z"+snapstr+".sovcirc"
sogrpname = fbase+modelname+"/so_z"+snapstr+".sogrp"

halos = genfromtxt(soname, names=True)
# gals = genfromtxt(galname)
# mstars = []
# for gal in gals: mstars.append(log10(gal[4] * UNIT_MASS * 1.e10 / HPARAM))
# mstars.append(-inf)
# mstars = array(mstars)
# hostmask = (halos['Msub'] > 0) & (mstars > 0)
# mstars = mstars[hostmask]
# msubs = log10(halos[hostmask]['Msub'] / HPARAM)
# # mfits = interp1d(msubs[::10], mstars[::10], kind="quadratic")
# mask = (msubs > 11.0) & (mstars > 8.0)
# msubs = msubs[mask]
# mstars = mstars[mask]
# coefs = polyfit(msubs, mstars, 2)
# x = linspace(10., 13.5, 100)
# y = polyval(coefs, x)

# plt.plot(msubs[::3], mstars[::3], "b.", alpha=0.4)
# plt.plot(x, y, "r-")


sohid = ioformat.rcol(sogrpname, [0], [0], linestart=1)

hf = h5py.File(snapname, "r")
print ("File Read.")

header = hf['Header'] # Read the header
Ngas = header.attrs['NumPart_Total'][0]
Ndark = header.attrs['NumPart_Total'][1]
Nstar = header.attrs['NumPart_Total'][4]
gp = hf['PartType0'] # Gas Particles
gp.keys()

print ("Reading HDF5: ", snapname)
# Get gas properties
mass = array(gp['Masses']) # Mass
u = array(gp['InternalEnergy']) # Temperature in K
ne = array(gp['ElectronAbundance'])
sfr = array(gp['StarFormationRate'])

sp = hf['PartType4']
massstar = array(sp['Masses'])

MeanWeight = (1.0 + 4.0 * XHE) / (1.0 + ne + XHE)
logT = log10(u * UNIT_T * MeanWeight)

mcold, mhot, mism = [0.0] * len(halos), [0.0] * len(halos), [0.0] * len(halos)
mstar = [0.0] * len(halos)

print ("Calculating for haloes: ")

for i in range(Ngas):
    hidx = sohid[i] - 1
    if(hidx < 0): continue
    if(sfr[i]): mism[hidx] += mass[i]
    else:
        if(logT[i] < TMAX): mcold[hidx] += mass[i]
        else: mhot[hidx] += mass[i]
for i in range(Nstar):
    offset = i + Ngas + Nstar
    hidx = sohid[offset] - 1
    if(hidx >= 0): mstar[hidx] += massstar[i]
    

mcold = log10(array(mcold) * 1.e10 / HPARAM)
mhot = log10(array(mhot) * 1.e10 / HPARAM)
mism = log10(array(mism) * 1.e10 / HPARAM)
mstar = log10(array(mstar) * 1.e10 / HPARAM)
        
print ("Writing result."        )
outname = fbase+modelname+"/so_z"+snapstr+".mbar"+suffix

print ("Output: ", outname)
fout = open(outname, "w")
fout.write("#Mvir Msub Mcold Mhot Mism Mstar\n")
for i in range(len(halos)-1):
    fout.write("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n" % \
               (log10(halos[i]['Mvir']/HPARAM), log10(halos[i]['Msub']/HPARAM),\
                mcold[i], mhot[i], mism[i], mstar[i]))
fout.close()
print ("DONE")
