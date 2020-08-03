from mymod import *
from cosmology import acosmic, tcosmic
import h5py
import os

# ================================
# Use showhalos.py to select haloes
# ================================

# Get info at each snapshot: logrho, logT, Halo Mvir
hparam = 0.7
omegab = 0.045
lbox = 25000.
UNIT_M = 1.989e43
UNIT_L = ac.kpc
UNIT_V = 1.e5
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
unit_length = lbox * ac.kpc
unit_time = sqrt(8.*pi/3.) * ac.mpc / (100.e5 * hparam)
unit_density = UNIT_M / (UNIT_L ** 3) / (1.8791e-29 * omegab)
unit_temp = UNIT_V ** 2 * (2./3.) * pc.mh / pc.k

model = "l25n288-phew-m5-spl"
haloid = 3357
snapi = 33
zstr = ("000"+str(snapi))[-3:]
datadir = "/proj/shuiyao/" + model + "/"
fbase = "/scratch/shuiyao/scidata/gadget3io/" + model + "/"
if(not os.path.exists(fbase + "haloparts/")):
    os.mkdir(fbase + "haloparts/")

HaloPos = array([0., 0., 0.])
HaloMass = 0.0
HaloRad = 0.0

def get_halo_info(hid):
    galname = datadir + "gal_z"+zstr+".stat"
    soname = datadir + "so_z"+zstr+".sovcirc"
    print "Galaxy Catalog: ", galname
    print "Halo Catalog: ", soname    
    x, y, z = ioformat.rcol(galname, [18, 19, 20])
    HaloPos[0] = (x[hid-1] + 0.5) * lbox
    HaloPos[1] = (y[hid-1] + 0.5) * lbox
    HaloPos[2] = (z[hid-1] + 0.5) * lbox
    msub, rsub = ioformat.rcol(soname, [6, 7], linestart=1)
    HaloMass = log10(msub[hid-1] / hparam)
    HaloRad = rsub[hid-1] / hparam
    print "Halo Information: "
    print "----------------"
    print "Pos = [%7.1f, %7.1f, %7.1f] (/h kpc)" % (HaloPos[0], HaloPos[1], HaloPos[2])
    print "Log(Mass/Msolar) = %6.3f" % (HaloMass)
    print "Radius = %6.2f kpc" % (HaloRad)
    return HaloPos, HaloMass, HaloRad

# pidset = set(ids)

snapname = datadir + "snapshot_"+zstr+".hdf5"
halostr = "h"+("00000" + str(haloid))[-5:]
outname = fbase + "haloparts/" + halostr + "_" + zstr
print "Reading: ", snapname
print "Writing: ", outname
HaloPos, HaloMass, HaloRad = get_halo_info(haloid)

hf = h5py.File(snapname, "r")
gp = hf['PartType0']
gp_pos = array(gp['Coordinates'])
# gp_rho = array(gp['Density'])
gp_temp = array(gp['InternalEnergy'])
gp_mass = array(gp['Masses'])
gp_ne = array(gp['ElectronAbundance'])
gp_mcloud = array(gp['PhEWMcloud'])
gp_hsml = array(gp['SmoothingLength'])
gp_sfr = array(gp['StarFormationRate'])

# Loop through all snapshots and find the particle
print "Start Searching Particles."
fout = open(outname, "w")
fout.write("%8.2f %8.2f %8.2f %6.3f %6.2f\n" % \
           (HaloPos[0], HaloPos[1], HaloPos[2], HaloMass, HaloRad))
fout.write("----------------\n")
fout.write("#x y z Hsml logT Mc SFflag\n")
ncount, nphew = 0, 0
for i in range(len(gp_pos)):
    if(not (i % (1000000))): print "Doing ", ("00000000" + str(i))[-8:]
    if(abs(gp_pos[i][0] - HaloPos[0]) > 1.2 * HaloRad): continue
    if(abs(gp_pos[i][1] - HaloPos[1]) > 1.2 * HaloRad): continue
    if(abs(gp_pos[i][2] - HaloPos[2]) > 1.2 * HaloRad): continue        
    # logrho = log10(gp_rho[i] * unit_density)
    meanweight = (1.0 + 4. * XHE) / (1. + XHE + gp_ne[i])
    logt = log10(gp_temp[i] * unit_temp * meanweight)
    hsml = gp_hsml[i]
    if(gp_sfr[i] > 0): sfflag = 1
    else: sfflag = 0
    mc = gp_mcloud[i]
    if(mc > 0): nphew += 1
    line = "%8.2f %8.2f %8.2f %7.2f %6.4f % 7.5f %d\n" % \
           (gp_pos[i][0], gp_pos[i][1], gp_pos[i][2], hsml, logt, mc, sfflag)
    fout.write(line)
    ncount += 1
fout.close()
print "Total Number of Gas Particles: ", ncount
print "Total Number of PhEW Particles: ", nphew
print "Done: ", outname
 
# print "DONE."
