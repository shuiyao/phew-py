# import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace, zeros
import h5py
import ioformat
import sys
# Based on readhdf5.py

# Calculate the fraction of winds in different components:
# Use l50n288.
# Halo Gas: Hot/Cold

errormsg = "Usage: calc_wind_fraction.py modelname lbox phewon(bool)"
if(len(sys.argv) != 4):
    raise ValueError(errormsg)
else:
    modelname = sys.argv[1]
    lbox = (float)(sys.argv[2])
    PhEWON = (int)(sys.argv[3])

redtab = genfromtxt("../redshifts.txt", names=True)

# BOXSIZE = 50000
# UNIT_MASS = 3469581.88
BOXSIZE = lbox * 1000.
UNIT_MASS = 3469581.88 * (lbox / 50.) ** 3
HPARAM = 0.7
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
UNIT_T = 1.e10 * (2./3.) * 1.6726e-24 / 1.3806e-16
MRES_HALO = 9.3e7 * 32

TTHRESH = 5.0

# UNIT_MASS = 3469581.88 * (12./50.) ** 3

fbase = "/nas/astro-th-nas/shuiyao/"
#list_of_snapshots = [58, 78, 98, 108]
list_of_snapshots = range(10, 108, 5)
#list_of_snapshots = [108]
fout = open("/home/shuiyao_umass_edu/scidata/"+modelname+"/"+modelname+".wfrac", "w")
fout.write("#z cold hot igm metalc metalh metali mcold mhot migm\n")

for snapnum in list_of_snapshots:
    snapstr = ("000"+str(snapnum))[-3:]
    snapname = fbase+modelname+"/snapshot_"+snapstr+".hdf5"
    sogrpname = fbase+modelname+"/so_z"+snapstr+".sogrp"
    soname = fbase+modelname+"/so_z"+snapstr+".sovcirc"
    zred = redtab[snapnum]['zred']

    halos = genfromtxt(soname, names=True)
    mvir = halos['Mvir'] / HPARAM
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
    metals = array(gp['Metallicity']).T[0] # Metals
    if(PhEWON):
        pkey = array(gp['PhEWKey'])
        wmass = array(gp['PhEWWindMass']) # WindMass
    else:
        pkey = zeros(len(mass))
    MeanWeight = (1.0 + 4.0 * XHE) / (1.0 + ne + XHE)
    logT = log10(u * UNIT_T * MeanWeight)

    print ("Calculating Wind Fractions.")

    nh = len(halos)
    m_igm, mw_igm, m_cold, mw_cold, m_hot, mw_hot = \
        0., 0., 0., 0., 0., 0.
    # Avoid divided by zero error
    mz_igm, mz_cold, mz_hot = 0., 0., 0.
    for i in range(Ngas):
        hidx = sohid[i] - 1
        if(sfr[i] > 0 and pkey[i] == 0):
            continue # ISM particles
        if(hidx < 0 or mvir[hidx] < MRES_HALO): # Non-halo gas
            m_igm += mass[i]
            mw_igm += wmass[i]
            mz_igm += mass[i] * metals[i]
        else:
            if(logT[i] < TTHRESH):
                m_cold += mass[i]
                mw_cold += wmass[i]
                mz_cold += mass[i] * metals[i]                
            else:
                m_hot += mass[i]
                mw_hot += wmass[i]
                mz_hot += mass[i] * metals[i]                

    mw_cold = mw_cold / m_cold if m_cold > 0 else 0.
    mw_hot = mw_hot / m_hot if m_hot > 0 else 0.
    mw_igm = mw_igm / m_igm if m_igm > 0 else 0.
    mz_cold = mz_cold / m_cold if m_cold > 0 else 0.
    mz_hot = mz_hot / m_hot if m_hot > 0 else 0.
    mz_igm = mz_igm / m_igm if m_igm > 0 else 0.
    line = "%5.3f %7.5e %7.5e %7.5e %7.5e %7.5e %7.5e %7.5e %7.5e %7.5e\n" % \
               (zred, mw_cold, mw_hot, mw_igm,
                mz_cold, mz_hot, mz_igm,
                m_cold, m_hot, m_igm)
    fout.write(line)
    print (line)

fout.close()
print ("Done.")

# outname = fbase+modelname+"/so_z108.mgas"
# fout = open(outname, "w")
# fout.write("#Mvir Msub Mcold Mhot Mism\n")
# for i in range(len(halos)-1):
#     fout.write("%6.3f %6.3f %6.3f %6.3f %6.3f\n" % \
#                (log10(halos[i]['Mvir']/HPARAM), log10(halos[i]['Msub']/HPARAM),\
#                 mcold[i], mhot[i], mism[i]))
# fout.close()
