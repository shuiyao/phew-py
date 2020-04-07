# Show some postshock quantities
# x300v1000, x300v1700, x1000v1700, x3000v3000, x300v3000, x300v1700c02, x300v1700c005

from mymod import *
from numpy import array
import conduction
import matplotlib.pyplot as plt
import config_mpl

rho_0 = array([2.7e-27, 2.7e-27, 8.1e-28, 2.7e-28, 2.7e-27, 2.7e-27, 2.7e-27])
T_0 = array([3.0e6, 3.0e6, 1.0e7, 3.0e7, 3.0e6, 3.0e6, 3.0e6])
rho_1 = array([4.4e-26, 4.6e-26, 1.5e-26, 8.4e-27, 9.7e-26, 3.5e-26, 2.9e-26])
T_1 = array([4.3e6, 1.2e7, 1.2e7, 2.7e7, 2.4e7, 1.4e7, 1.9e7])
machs = array([3.8, 6.5, 3.5, 3.6, 11.4, 6.5, 6.5])
clrs = ['blue', 'red', 'green', 'cyan', 'magenta', 'orange', 'yellow']
models = [r"$\chi300v1000$", r"$\chi300v1700$", r"$\chi1000v1700$", r"$\chi3000v3000$", r"$\chi300v3000$", r"$\chi300v1700c5$", r"$\chi300v1700c20$"]
sizes = [150, 150, 150, 150, 150, 150]

ratio_rho = rho_1 / rho_0
ratio_T = T_1 / T_0

fig = plt.figure(1, figsize=(8,7))
ax = fig.add_subplot(111)
syms = []
for i in range(len(ratio_rho)):
    sym = ax.scatter(ratio_rho[i], ratio_T[i], marker="*", s=sizes, color=clrs[i])
    syms.append(sym)

# Mode 1: Mach from simulation, q range from 0 to 1
qs = linspace(0.01, 0.99, 100)
linei = 0
for mach in machs:
    fac_rho, fac_T = [], []
    for q in qs:
        fac_rho.append(conduction.ratio_density(q, mach))
        fac_T.append(conduction.ratio_temperature(q, mach))
    plt.plot(fac_rho, fac_T, "--", color=clrs[linei])
    linei += 1

# qs = linspace(0.01, 0.99, 10)
machs = linspace(1.0, 11.0, 6)
for mach in machs:
    fac_rho, fac_T = [], []    
    qs = linspace(0.01, 0.99, 100)    
    for q in qs:
        fac_rho.append(conduction.ratio_density(q, mach))
        fac_T.append(conduction.ratio_temperature(q, mach))
    plt.plot(fac_rho, fac_T, "k:")
qs = array([0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0])
txtoffset = 0.0
firsttxt = True
for q in qs:
    if(firsttxt):
        plt.text(fac_rho[0]-1.2, fac_T[0], "$\hat{q}_s = $", fontsize=14, color="black")
        firsttxt = False
    if q == 1.0: q = 0.99
    fac_rho, fac_T = [], []
    machs = linspace(1.0, 11.0, 100)    
    for mach in machs:
        fac_rho.append(conduction.ratio_density(q, mach))
        fac_T.append(conduction.ratio_temperature(q, mach))
    plt.plot(fac_rho, fac_T, "k:")
    if(q == 0.99): 
        plt.text(fac_rho[-1]+txtoffset, fac_T[-1], str(q)[:4], fontsize=14, color="black")
    else:
        plt.text(fac_rho[-1]+txtoffset, fac_T[-1], str(q)[:3], fontsize=14, color="black")        

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\rho_{ps}/\rho_a$")
plt.ylabel(r"$T_{ps}/T_a$")
plt.legend(syms, models, loc='upper right')
plt.savefig("qs.pdf")
plt.show()
