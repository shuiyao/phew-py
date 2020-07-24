# Find the fion from tabulated data.

from numpy import genfromtxt, sqrt, log10
from astroconst import pc, ac

# Ion lookup table definitions
MAXIONS = 9 # 35
NHPTS = 240 # Number of n_h gridpoints in lookup table
TPTS = 140 # Number of T gridpoints in lookup table
GALPTS = 11
# Table limits
NHLOW = -9.0
DELTANH = 0.05
TLOW = 2.5
DELTAT = 0.05
GALLOW = 6.0
DELTAGAL = 0.5

ionfolder = "/scratch/shuiyao/specexbin/ionfiles/"

Ions = genfromtxt(ionfolder+"specions_i9.dat", dtype=("S8,f8,f8,f8,f8,i8,f8"),
                  names = ("name","wvlen","Xsec","atomwt","fraction","zcol","alpha"))

Ions['Xsec'] *= 2.648e-2 * Ions['wvlen'] * 1.e-13
bsys = sqrt(2. * pc.k / (pc.mh * Ions['atomwt'])) / 1.e5

tab = genfromtxt(ionfolder+"lt02HM12_i9", dtype=("f8,f8,f8,f8,f8,f8,f8,f8,f8"),
                 names = tuple(Ions['name']))
tab = tab.reshape(NHPTS, TPTS)
# Now use ---- tab[inh][itemp][ion] ---- to get the value

nhl = [0.0] * NHPTS
tl = [0.0] * TPTS
# fractab = [[[0.0] * TPTS] * NHPTS] * 9
for inh in range(NHPTS): nhl[inh] = NHLOW + inh * DELTANH
for itemp in range(TPTS): tl[itemp] = TLOW + itemp * DELTAT

def IonFrac(temp, density, ionid):
    nhmax, nhmin, tmax, tmin = 0.000001, 999, 11., 1.e7

    n_h = density / pc.mh
    if(n_h<nhmin): nhmin=n_h
    if(n_h>nhmax): nhmax=n_h
    if(temp<tmin): tmin=temp
    if(temp>tmax): tmax=temp

    # n_h is cgs number density of Hydrogen
    # temp is temperature in K
  
    inh = (int)((log10(n_h) - NHLOW)/DELTANH)
    itemp = (int)((log10(temp) - TLOW)/DELTAT)

    if(inh>NHPTS-2): inh=NHPTS-2
    if(itemp>TPTS-2): itemp=TPTS-2
    if(inh<0): inh=0
    if(itemp<0): itemp=0
    t = (log10(n_h)-nhl[inh])/DELTANH
    u = (log10(temp)-tl[itemp])/DELTAT

    y1 = tab[inh][itemp][ionid]
    y2 = tab[inh+1][itemp][ionid]
    y3 = tab[inh+1][itemp+1][ionid]
    y4 = tab[inh][itemp+1][ionid]
    ionfraction = 10. ** ((1.-t)*(1.-u)*y1+t*(1.-u)*y2+t*u*y3+(1.-t)*u*y4)
    return ionfraction

print "DONE"
