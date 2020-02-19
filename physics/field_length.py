# Things that relate to the Field Length

import bowshock
import conduction
from astroconst import pc, ac
from numpy import sqrt, logspace, linspace
from numpy import genfromtxt

conversion_T_keV = pc.k / pc.ev / 1.e3
HIz0, HIz2 = 8.89e-14, 3.74e-12
HIz0 *= pc.ev
HIz2 *= pc.ev

# Sharma et al. 2010, Tozzi & Norman 2001
def lambda_cooling(T):
    T = T * conversion_T_keV
    if(T < 0.0017235):
        return 1.544e-22 * (T / 0.0017235) ** 6
    elif(T < 0.02):
        return 6.72e-22 * (T / 0.02) ** 0.6

def tau_cool(n, T, gamma=5./3.):
    return 1. / (gamma - 1.) * pc.k * T / (n * lambda_cooling(T))

def mfp(n, T): # McKee & Begelman 1990
    T7 = T / 1.e7
    return 0.284 * T7 * T7 / n * ac.pc

def kappa_sat(nh, Th): # McKee & Begelman 1990
    y = 1.19 * nh * pc.k * mfp(nh, Th) * \
        sqrt(pc.k * Th / pc.me)
    return y

def l_field(nc, Tc, nh, Th, saturation=False):
    if(saturation == False):
        y = conduction.kappa_cond(nh, Th) * Th / (nc * nc * lambda_cooling(Tc))
    else:
        y = kappa_sat(nh, Th) * Th / (nc * nc * lambda_cooling(Tc))
    return sqrt(y)

def l_field_eq(nc, nh, Th, H_HI=HIz0): # Use equilibrium T, that is, heating rate
    y = conduction.kappa_cond(nh, Th) * Th / (nc * H_HI)
    return sqrt(y)

nc, Tc = 1.0, 1.e4
nh, Th = 1.e-3, 3.e6

print "L_f = %5.3f [pc]" % (l_field(nc, Tc, nh, Th) / ac.pc)
print "L_f(sat) = %5.3f [pc]" % (l_field(nc, Tc, nh, Th, saturation=True) / ac.pc)
print "L_mfp = %5.3f [pc]" % (mfp(nh, Th) / ac.pc)
tab = genfromtxt("hm12.dat", usecols=(0, 2), names=['z', 'HI_h'])

import matplotlib.pyplot as plt
def show_cooling_heating(nc=1.0):
    Ts = linspace(1.e3, 2.e4, 100)
    ratio = []
    for Tc in Ts:
        ratio.append(lambda_cooling(Tc) * nc / HIz0)
    plt.plot(Ts, ratio, "k-")
    ratio = []
    for Tc in Ts:
        ratio.append(lambda_cooling(Tc) * nc / HIz2)
    plt.plot(Ts, ratio, "k--")
    plt.axhline(1.0, color="blue")
    plt.yscale("log")
    plt.show()
    
nc, Tc = 10.0, 1.e4
nh = 1.e-3
Thot = logspace(5.5, 7.5, 100)
Lfs, Lfeq, Lmfps = [], [], []
for Th in Thot:
    Lfs.append(l_field(nc, Tc, nh, Th) / ac.pc)
    Lfeq.append(l_field_eq(nc, nh, Th, HIz2) / ac.pc)    
    Lmfps.append(mfp(nh, Th) / ac.pc)
plt.plot(Thot, Lfs, "b-")
plt.plot(Thot, Lfeq, "g-")
plt.plot(Thot, Lmfps, "r-")
plt.axhline(y=20.0, color="black", linestyle=":")
plt.xscale("log")
plt.yscale("log")
plt.show()

