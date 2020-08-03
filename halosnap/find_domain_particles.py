from mymod import *
from cosmology import acosmic, tcosmic
import h5py
import os

# ================================
# Use showhalos.py to select haloes
# ================================

model = "l12n144-phew-movie"
lbox = 12000.
haloid = 2318 # log(Msub) ~ 12.0 Spiral

# haloid = 1503 # log(Msub) ~ 11.0
#haloid = 1972 # log(Msub) ~ 12.0
#haloid = 251 # log(Msub) ~ 13.0

# Get info at each snapshot: logrho, logT, Halo Mvir
hparam = 0.7
omegab = 0.045
UNIT_M = 1.989e43
UNIT_L = ac.kpc
UNIT_V = 1.e5
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
unit_length = lbox * ac.kpc
unit_time = sqrt(8.*pi/3.) * ac.mpc / (100.e5 * hparam)
unit_density = UNIT_M / (UNIT_L ** 3) / (1.8791e-29 * omegab)
unit_temp = UNIT_V ** 2 * (2./3.) * pc.mh / pc.k

fbase = "/scratch/shuiyao/scidata/halosnap/" + model + "/"
galname = "/proj/shuiyao/" + model + "/gal_z200.stat"
x, y, z = ioformat.rcol(galname, [18, 19, 20])
if(not os.path.exists(fbase + "box/")):
    os.mkdir(fbase + "box/")

# HaloPos = array([0.0, 0.1, 0.12]) # The center of the domain
# HaloPos = array([0.14268, 0.12580, -0.16811]) # 1503
HaloPos = [x[haloid-1], y[haloid-1], z[haloid-1]]
HaloPos = array(HaloPos)
print "Galaxy Position: ", HaloPos
HaloPos = (HaloPos + 0.5) * lbox

HaloRad = 0.025 # The (half) box size
HaloRad *= lbox

halostr = ("00000"+str(haloid))[-5:]
datadir = "/proj/shuiyao/" + model + "/"

# Loop Begins Here
# ================================

for snapi in range(201):
    if(snapi < 10): continue
    zstr = ("000"+str(snapi))[-3:]
    snapname = datadir + "snapshot_"+zstr+".hdf5"

    outname = fbase + "movie/" + "box_" + halostr + "_" + zstr
    print "Reading: ", snapname
    print "Writing: ", outname

    fout = open(outname, "w")
    fout.write("%8.2f %8.2f %8.2f %6.3f %6.2f\n" % \
               (HaloPos[0], HaloPos[1], HaloPos[2], 0.0, HaloRad))
    fout.write("----------------\n")

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
    fout.write("#x y z Hsml logT Mc SFflag\n")
    ncount, nphew = 0, 0
    for i in range(len(gp_pos)):
        if(not (i % (1000000))): print "Doing ", ("00000000" + str(i))[-8:]
        if(abs(gp_pos[i][0] - HaloPos[0]) > 1.0 * HaloRad): continue
        if(abs(gp_pos[i][1] - HaloPos[1]) > 1.0 * HaloRad): continue
        if(abs(gp_pos[i][2] - HaloPos[2]) > 1.0 * HaloRad): continue        
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
 
print "DONE."
