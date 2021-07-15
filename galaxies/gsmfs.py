# OPTIMIZED FOR SMF

import ioformat
import matplotlib.pyplot as plt
from numpy import histogram, logspace, linspace, array, sqrt, log10
from myinit import *

# Compare SMF at different redshift.
# Write smf_###.txt file

DATABASEDIR = DIRS['DATA']
# modelname = "l50n576-phew-m5"
modelname = "l50n288-phewoff-sph"
SKIDBASE = modelname+"/"
#SKIDBASE = "p25n144gwl_z1b/"

unit_l = 50000.
unit_m = 3469578.81574
unit_v = 1727.47074736

class groups():
    def __init__(self, filename, t=1., unit_system="standard"):
        self.filename = filename
        self.L = 1.
        self.t = t
        if unit_system=="gadget":
            self.unit_system = "gadget"
            self.unit_l = unit_l
            self.unit_m = unit_m
            self.unit_v = unit_v
        elif unit_system=="tipsy":
            self.unit_system = "tipsy"
            self.unit_l = 1.
            self.unit_m = 1.
            self.unit_v = 1.
        elif unit_system=="standard":
            self.unit_system = "standard"
            self.unit_l = unit_l
            self.unit_m = unit_m * 1.e10
            self.unit_v = unit_v * t#sqrt(self.t)

        f = open(filename, "r")
        cols = ioformat.readcol(f, 2)
        f.close()
        self.L = self.L * self.unit_l * 0.001 #Mpc
        self.gid = cols.cols[0]
        self.n = cols.cols[1]
        self.masstot = array(cols.cols[2]) * self.unit_m
        self.massgas = array(cols.cols[3]) * self.unit_m
        # self.massstar = array(cols.cols[4]) * self.unit_m * 7.68 / 1.23
        self.massstar = array(cols.cols[4]) * self.unit_m
        self.vcirmax = array(cols.cols[5]) * self.unit_v
        self.vcirhalf = array(cols.cols[6]) * self.unit_v
        self.vcirouter = array(cols.cols[7]) * self.unit_v
        self.rvcirmax = array(cols.cols[8]) * self.unit_l
        self.rhalfmass = array(cols.cols[9]) * self.unit_l
        self.router = array(cols.cols[10]) * self.unit_l
        self.vdisp = array(cols.cols[11]) * self.unit_v

def read(filename, t=1.):
    return groups(filename, t)

def smf(infile, outfile, Nbins=40, Mmin=8., Mmax=13., mcorr=1.0, plot=False):
    grp = read(infile)
    bins = logspace(Mmin, Mmax, Nbins+1)
#   dlogM = dM / M
    dlogM = (Mmax - Mmin)/(Nbins+1)
    grp.massstar *= mcorr
    hist, edges = histogram(grp.massstar, bins=bins)#, normed=True)
    f = open(outfile, "w")
    f.write("# Header: Nbins, MassMin, MassMax\n")
    f.write(str(Nbins)+" "+str(10.**Mmin)+" "+str(10.**Mmax)+"\n")
    mid, phi = [], []
    for i in range(Nbins):
        mid.append(10.**((log10(edges[i])+log10(edges[i+1]))/2.))
#        dlogM = (edges[i+1] - edges[i]) / edges[i]
        phi.append(hist[i]/(grp.L)**3/dlogM)
        f.write(str(mid[i])+" "+str(phi[i])+"\n")
    f.close()
    if plot==True: 
        plt.plot(mid, phi, ".-")
        plt.xscale("log")
        plt.yscale("log")

def batch():
    # CHANGE SKIDBASE and outputfile
    # for n in ["033","058","043"]:
    for n in [snapstr]:
        infile = DATABASEDIR+SKIDBASE+"gal_z"+n+".stat"
        # outfile = "./smf_p50n576gw_"+n+".txt"
        # outfile = "./smf_p12n144gw_"+n+".txt"
        # outfile = "./smf_r32n512_"+n+".txt"
        print (infile+" >")
        print (outfile)
        smf(infile, outfile, plot=True)

for snapstr in ["033", "058", "078", "108", "098"]:    
    outfile = DIRS['SCIDATA']+modelname+"/gsmf_"+snapstr+".txt"o
    batch()

