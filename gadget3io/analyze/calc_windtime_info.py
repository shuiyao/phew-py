from mymod import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic

# ------------ Purpose ------------
# Write .windage data
# Each row: the amount of wind mixed at time a
# a, mh11, mh12, mh13, all

hparam = 0.7
unit_m = 1.e10 / hparam

CUMULATIVE_DISTRIBUTION = False

nbins = 100
aedges = linspace(0., 1., nbins+1)
acen = 0.5 * (aedges[1:] + aedges[:-1])
da = 1. / nbins

def find_windtime_bin(awind):
    bin_idx = awind / da
    return (int)(bin_idx)

models = ["l25n144-phew-m5-spl", "l25n288-phew-m5-spl"]

REDSHIFT = 0.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"
if(REDSHIFT == 0.25): zstr = "098"
ALIM = acosmic(tcosmic(1./(1.+REDSHIFT)) - 1.e9)

mmin = [11.0, 11.85, 12.85, 0.0]
mmax = [11.5, 12.15, 13.15, 16.0]

for modeli in range(len(models)):
    fname = "/home/shuiyao_umass_edu/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".starinfo"
    soname = "/nas/astro-th-nas/shuiyao/"+models[modeli]+"/so_z"+zstr+".sovcirc"
    foutname = "/home/shuiyao_umass_edu/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".windage"
    print ("Doing: ", fname)
    stars = genfromtxt(fname, names=True)
    halos = genfromtxt(soname, names=True)
    halos['Msub'] = log10(halos['Msub'] / 0.7)
    nstars = len(stars)
    stars = stars[stars['WindMass'] > 0.0]
    # stars = stars[stars['a_acc'] > 0.80] # z < 0.25
    stars = stars[stars['a_acc'] > ALIM] # Winds accreted within 1 Gyr
    print ("%d out of %d stars selected. " % (len(stars), nstars))
    windage = []
    for mi in range(4):
        mbins_totw = array([0.0] * nbins)
        for i, s in enumerate(stars):
            hidx = (int)(s['HID'] - 1)
            if(mmin[mi] < halos[hidx]['Msub'] < mmax[mi]):
                aidx = find_windtime_bin(s['WindAge'])
                mbins_totw[aidx] += s['WindMass']
        windage.append(mbins_totw)
    fout = open(foutname, "w")
    fout.write("#atime Mh11 Mh12 Mh13 All\n")
    for i in range(nbins):
        line = "%7.5f %7.5e %7.5e %7.5e %7.5e\n" % \
            (acen[i], windage[0][i], windage[1][i], windage[2][i], windage[3][i])
        fout.write(line)
    fout.close()
        
print ("Done.")
