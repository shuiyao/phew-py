import phew
from astroconst import pc, ac
from numpy import log10, log, sqrt, array
import matplotlib.pyplot as plt

GAMMA = 5./3.
MCLOUD = 2.0e38

import config_mpl

def remove_spurious_particles(track):
    dr = track['dr']
    atime = track['atime']
    vrel = track['vrel']
    # if(max(dr) >= 1000.): return 0
    if(dr[0] >= 10.): return 0
    for i in range(len(dr)-1):
        if(atime[i+1] - atime[i] > 0.01): return 0
        if(abs(dr[i+1] - dr[i]) > 0.1 * dr[-1]): return 0
        if(abs(vrel[i+1] - vrel[i]) > 0.1 * vrel[0]): return 0
    return 1

def draw(PhEWParticles, ax, fieldc, \
         nskip=1, color_selection=True, \
         color_min=None, color_max=None, logcscale=False, alpha=0.2):
    if(color_min == None): color_min = min(fieldc)
    if(color_max == None): color_max = max(fieldc)
    phew.set_colors_phew_particles(PhEWParticles, fieldc, color_min, color_max, logscale=logcscale)
    iskip = 0
    for PhEWP in PhEWParticles:
        if(iskip < nskip):
            iskip += 1
            continue
#        flag = remove_spurious_particles(PhEWP.track)
        flag = 1
        if(flag == 0): continue
        if(color_selection == True):
            if(not (color_min < PhEWP.track[fieldc][0] < color_max)):
                continue
        cs_a = sqrt(GAMMA * pc.k * PhEWP.track['T_a'] / (0.60 * pc.mh))
        a0 = PhEWP.track['atime'][0]
        # xarr = PhEWP.track['vrel'] / cs_a
        # yarr = PhEWP.track['T_c'] / PhEWP.track['T_a']
        xarr = PhEWP.track['Mach']
        yarr = PhEWP.track['M_c'] / MCLOUD
        # xarr = PhEWP.track['atime'] - a0
        # yarr = PhEWP.track['dr']
        ax.plot(xarr, yarr, "-", color=PhEWP.color, alpha=alpha)
        ax.plot(xarr[0], yarr[0], "^", color=PhEWP.color, markersize=8)
        ax.plot(xarr[-1], yarr[-1], "*", color=PhEWP.color, markersize=8)
        iskip = 0

filename = "/proj/shuiyao/m6n64beta6/WINDS/z2/sorted.phews"
#filename = "/proj/shuiyao/l25n144phew/WINDS/z2/sorted.phews"
Fields = phew.PhEWFields()
info_fields = Fields.get_field_info(["atime","ID","dr","dv","vrel","Mach","M_c","rho_a","Mstar","T_c","T_a"], verbose=True)
tab = phew.read_phew_fields(filename, info_fields)
# subtab = tab[tab['T_c'] < 1.e4]
pparts = phew.create_phew_particles(tab)

fig = plt.figure(1, figsize=(6,6))
ax = fig.add_subplot(111)
draw(pparts, ax, "T_a", nskip=10, logcscale=True, color_selection=True, color_min=2.e4, color_max=1.e7, alpha=0.3)
ax.set_xlim(0.0, 10.0)
ax.set_ylim(0.0, 1.0)
ax.plot([1.0, 1.0], [0.0, 1.0], "k:") # vertical line
ax.plot([0.0, 10.0], [0.05, 0.05], "k:")
# ax.set_ylim(1.e-3, 2.0)
# ax.set_yscale("log")
ax.set_xlabel(r"$\mathcal{M}$")
ax.set_ylabel(r"$M_c/M_c(0)$")
plt.show()
