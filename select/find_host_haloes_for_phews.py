import matplotlib.pyplot as plt
import h5py
import ioformat
import sys
from numpy import genfromtxt, log10

# Find the host haloes for PhEW particles.
# Input: initwinds.*, .sogrp, particleIDs and masses
# Output: idx PhEWKey HID Msub

errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) modelname"
if(len(sys.argv) != 4):
    print errormsg
    sys.exit(1)
else:
    NCPU = int(sys.argv[1])
    REDSHIFT = float(sys.argv[2])
    MODELNAME = sys.argv[3]

# NCPU = 128
# MODELNAME = "l25n144-phew"
# REDSHIFT = 1.0

folder = "/proj/shuiyao/"
if(REDSHIFT == 0.2):
    snapnum = 100
    zstr = 'z0'
if(REDSHIFT == 1.0):
    snapnum = 78
    zstr = 'z1'
if(REDSHIFT == 2.0):
    snapnum = 58
    zstr = 'z2'
if(REDSHIFT == 4.0):
    snapnum = 33
    zstr = 'z4'
    
snapstr = ("000"+str(snapnum))[-3:]
fbase = folder + MODELNAME
snapname = fbase + "/snapshot_" + snapstr + ".hdf5"
sogrpname = fbase + "/so_z" + snapstr + ".sogrp"
sovcircname = fbase + "/so_z" + snapstr + ".sovcirc"
# grpname = fbase + "/gal_z" + snapstr + ".grp"
outname = fbase + "/snapshot_" + snapstr + ".phews"

def number_of_lines(filename):
    f = open(filename, "r")
    count = 0
    for line in f:
        count += 1
        if(count >= 2): return 2
    return count

def load_files():
    global pids
    global mass
    global hid
    global attrs
    global haloes
    
    print "Loading Files... "
    hf = h5py.File(snapname, "r")

    header = hf['Header']
    attrs = header.attrs
    attrs.keys()
    gp = hf['PartType0']
    gp.keys()
    pids = gp['ParticleIDs']
    mass = gp['Masses']

    hid = ioformat.rcol(sogrpname, [0], [0])
    haloes = genfromtxt(sovcircname, usecols=(1,2,6,7),
                        names=['Mvir','Rvir','Msub', 'Rsub'], skip_header=1)

# Construct hash map from initwinds.*
def build_maps():
    print "Building maps PID -> Mass, PhEWKey ... "
    phewbase = fbase + "/WINDS/" + zstr + "/"
    for icpu in range(NCPU):
        finitwinds = phewbase + "initwinds." + str(icpu)
        flag = number_of_lines(finitwinds)
        if(flag == 0): continue
        tab = genfromtxt(finitwinds, usecols=(1,2,17), dtype='i8, f8, i8')
        if(flag == 2): # More than two lines
            for line in tab:
                pid_to_key[line[2]] = line[0]
                pid_to_mass[line[2]] = line[1]
        else: # Only one line
            line = tab.item()
            pid_to_key[line[2]] = line[0]
            pid_to_mass[line[2]] = line[1]

# Match snapfiles to initwinds.*
pid_to_idx = dict()
def match_phews():
    print "Matching PhEWs ... "
    Nsph = attrs['NumPart_Total'][0]
    for i in range(Nsph):
        if(i % 100000 == 0): print "i = ", i
        pid = pids[i] # pids = gp['ParticleIDs']
        if(pid in pid_to_key): # PID found in initwinds.*
            phew_key = pid_to_key[pid]
            phew_mass = pid_to_mass[pid]
            if(pid in pid_to_idx):
                idx = pid_to_idx[pid]
                if(abs(mass[i] - phew_mass) < abs(mass[idx] - phew_mass)): # a better match
                    pid_to_idx[pid] = i
            else:
                pid_to_idx[pid] = i

# Write:
def write_file():
    fout = open(outname, "w")
    fout.write("#idx PhEWKey HID LogMvir LogMsub Rvir\n")
    for pid in pid_to_key:
        idx = pid_to_idx[pid]
        if(abs(hid[idx]) < len(haloes)): # Ignore last line
            halo = haloes[abs(hid[idx])-1]
            line = "%8d %10d %5d %6.3f %6.3f %6.1f\n" % (
                idx, pid_to_key[pid], hid[idx],
                log10(halo['Mvir']/0.7), log10(halo['Msub']/0.7), halo['Rvir']/0.7)
        fout.write(line)
    fout.close()


load_files()
pid_to_key = dict()
pid_to_mass = dict()
build_maps()
match_phews()
write_file()
