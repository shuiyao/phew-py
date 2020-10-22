# import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace
import h5py
import ioformat

# Based on readhdf5.py

# Calculate the fraction of winds from SFRINFO accretion info:
# Use l50n288.
# Halo Gas: Hot/Cold

# 1. Get stars.AccKey from z108
# 2. Find the accretion events (Should do a separate table)
# 3. Put accretion into redshift bins.

BOXSIZE = 50000
UNIT_MASS = 3469581.88
HPARAM = 0.7
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
UNIT_T = 1.e10 * (2./3.) * 1.6726e-24 / 1.3806e-16
MRES_HALO = 9.3e7 * 32

# UNIT_MASS = 3469581.88 * (12./50.) ** 3

# modelname = "l50n288-phewoff"
# modelname = "l50n288-phew-m5"
modelname = "p50n288fiducial"
fbase = "/nas/astro-th/shuiyao/"
snapname = fbase+modelname+"/snapshot_108.hdf5"
sogrpname = fbase+modelname+"/so_z108.sogrp"
soname = fbase+modelname+"/so_z108.sovcirc"

halos = genfromtxt(soname, names=True)
sohid = ioformat.rcol(sogrpname, [0], [0], linestart=1)
hf = h5py.File(snapname, "r")
print ("File Read: ", snapname)

header = hf['Header'] # Read the header
Ngas = header.attrs['NumPart_Total'][0]
gp = hf['PartType0'] # Gas Particles
gp.keys()

print ("Reading HDF5: ", snapname)
# Get gas properties
mass = array(gp['Masses']) # Mass
u = array(gp['InternalEnergy']) # Temperature in K
ne = array(gp['ElectronAbundance'])
sfr = array(gp['StarFormationRate'])
wmass = array(gp['WindMass']) # WindMass
MeanWeight = (1.0 + 4.0 * XHE) / (1.0 + ne + XHE)
logT = log10(u * UNIT_T * MeanWeight)

for i in range(Ngas):
    hidx = sohid[i] - 1
    if(sfr[i] > 0): continue
    if(hidx < 0 or halos[hidx]['Mvir']/HPARAM < MRES_HALO): # Non-halo gas
        m_igm += mass[i]
        mw_igm += wmass[i]
    else:
        if(logT[i] < 5.5):
            m_cold[hidx] += mass[i]
            mw_cold[hidx] += wmass[i]
        else:
            m_hot[hidx] += mass[i]
            mw_hot[hidx] += mwass[i]

print ("Writing result."        )
fout.write("%5.3f %7.5e %7.5e %7.5e\n" % \
           (zred, mw_cold/m_cold, mw_hot/m_hot, mw_igm/m_igm)

# outname = fbase+modelname+"/so_z108.mgas"
# fout = open(outname, "w")
# fout.write("#Mvir Msub Mcold Mhot Mism\n")
# for i in range(len(halos)-1):
#     fout.write("%6.3f %6.3f %6.3f %6.3f %6.3f\n" % \
#                (log10(halos[i]['Mvir']/HPARAM), log10(halos[i]['Msub']/HPARAM),\
#                 mcold[i], mhot[i], mism[i]))
# fout.close()
print ("DONE")
