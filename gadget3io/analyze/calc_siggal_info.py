from myinit import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic


hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 25
SIGMAX = 150.0
sigedges = linspace(0., SIGMAX, nbins+1)
sigcen = 0.5 * (sigedges[1:] + sigedges[:-1])
dsig = SIGMAX / nbins

ZSOLAR = log10(0.0122)

def find_siggal_bin(sig):
    if(sig >= SIGMAX): sig = SIGMAX - 1.0
    bin_idx = sig / dsig
    return (int)(bin_idx)

models = ["l25n144-phew-m5", "l25n288-phew-m5"]

REDSHIFT = 0.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"
if(REDSHIFT == 0.25): zstr = "098"
# ALIM = acosmic(tcosmic(1./(1.+REDSHIFT)) - 1.e9)
ALIM = 0.0

mmin = [11.0, 11.85, 12.85]
mmax = [11.5, 12.15, 13.15]

for modeli in range(len(models)):
    fname = DIRS['SCIDATA']+models[modeli]+"/"+models[modeli]+"_"+zstr+".starinfo"
    soname = DIRS['DATA']+models[modeli]+"/so_z"+zstr+".sovcirc"
    print "Doing: ", fname
    stars = genfromtxt(fname, names=True)
    halos = genfromtxt(soname, names=True)
    halos['Msub'] = log10(halos['Msub'] / 0.7)
    nstars = len(stars)
    stars = stars[stars['WindMass'] > 0.0]
    stars = stars[stars['a_acc'] > ALIM] # z < 0.25    
    print "%d out of %d stars selected. " % (len(stars), nstars)
    for mi in range(3):
        mbins_siggal = array([0.0] * nbins)
        for i, s in enumerate(stars):
            hidx = (int)(s['HID'] - 1)
            if(mmin[mi] < halos[hidx]['Msub'] < mmax[mi]):
                aidx = find_siggal_bin(s['WindSig'])
                mbins_siggal[aidx] += s['WindMass']
        
print ("Done.")
