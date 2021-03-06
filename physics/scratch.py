from astroconst import pc, ac
from numpy import exp, log, pi, linspace
from scipy.optimize import bisect
import matplotlib.pyplot as plt

def func(y, x0=1.0, t0=0.0):
    return exp(x0 * (y - 1.0)) * (1.0 + y * t) - 1.0

rho0 = 10. * pc.mh * 0.6
r_c = 20.*ac.pc
mlra = 4.46e-15 * 1.e7** 2 / (rho0 * r_c ** 2)
print log(mlra * 50. * ac.myr)

# https://www.dropbox.com/home/Ion%20Analysis
# Very crudely read from J'Neil's absorption images
# x300v1000:
# OVI 13 - 14;
# HI  13 - 18
# CIV 15 - 17
# x300v1700:
# HI  17 - 18
# CIV 15 - 16
# OVI 13 - 14

n_c = 1.0
rho_c = n_c * pc.mh * 0.62
r_c = 100 * ac.pc
M_c = 4./3. * pi * r_c ** 3 * rho_c
print "Mc = ", M_c / ac.msolar
# On average:
N_c = 4./3. * n_c * pi * r_c
print "N_c (avg.) = ", N_c
# Scale: n_c ~ r_c ^ (-3), N_c ~ n_c * r_c ~ r_c ^ (-2)
