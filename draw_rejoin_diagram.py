import phew
from astroconst import pc, ac
from numpy import log10, log, sqrt, array, isnan
import matplotlib.pyplot as plt
from numpy import genfromtxt
import matplotlib as mpl
from random import random

# Do match_initwinds_rejoin.py first

GAMMA = 5./3.
MCLOUD = 2.0e38
CMAP = "jet"
NSKIP = 1

import config_mpl

# model = "l25n144-phew-m5kh100fs10"
# model = "l25n144-phew-norecouple"
model = "l50n288-phew-m5-rec"
filename = "/nas/astro-th/shuiyao/"+model+"/WINDS/z1/sorted.phews"
fphewsname = "/home/shuiyao_umass_edu/scidata/"+model+"/phewsinfo.z1"

def select_particles(PhEWParticles, ntot=60, mmin=11.0, mmax=13.5):
    selected = []
    nbins = (int)(ntot / 2)
    mbins = [0] * nbins
    dm = (mmax - mmin) / (float)(nbins)
    for i, PhEWP in enumerate(PhEWParticles):
        idx = (PhEWP.mvir - mmin) / dm
        if(isnan(idx) or idx < 0.0 or idx >= nbins): continue
        idx = (int)(idx)
        if(mbins[idx] >= 2): continue
        mbins[idx] += 1
        selected.append(PhEWP)
    return selected

def select_particles_phew(PhEWParticles):
    selected = []
    random_norm = 0.5
    for i, PhEWP in enumerate(PhEWParticles):    
        random_cut = random_norm + 0.3 * (PhEWP.mvir - 12.0)
        if(random() > random_cut or PhEWP.mvir >= 15.0):
            continue
        selected.append(PhEWP)
    return selected

counter_spurious_particles = [0, 0, 0, 0, 0]
def remove_spurious_particles(track):
    dr = track['dr']
    atime = track['atime']
    vrel = track['vrel']
    # if(max(dr) >= 1000.): return 0
    if(dr[0] >= 50.): return 1
    for i in range(len(dr)-1):
        if(atime[i+1] - atime[i] > 0.01): return 2
        if(abs(dr[i+1] - dr[i]) > 0.2 * dr[-1]): return 3
        if(abs(vrel[i+1] - vrel[i]) > 0.2 * vrel[0]): return 4
    return 0
# One Example: [193, 20, 2, 215, 99]

def get_initwinds_and_rejoin_info(PhEWParticles, filename):
    # Note: some particles has HID = 0;
    # Their halo properties are from the last line of sovcirc. (Rvir < 0)
    tab = genfromtxt(filename, names=True, dtype=('f8,f8,f8,f8,f8,f8,f8,i8,f8,f8,f8,i8'))
    key_to_idx = dict()
    for i in range(len(tab)):
        if(tab[i]['Rvir'] > 0): key_to_idx[tab[i]['PhEWKey']] = i
    for PhEWP in PhEWParticles:
        if(PhEWP.key in key_to_idx):
            idx = key_to_idx[PhEWP.key]
            cs_a = sqrt(GAMMA * pc.k * 10.**tab[idx]['T_a'] / (0.60 * pc.mh))
            PhEWP.atimei = tab[idx]['a_i']
            PhEWP.atimer = tab[idx]['a_rejoin']
            PhEWP.machi = tab[idx]['Vinit'] * 1.e5 / cs_a # not working any more, the Vinit is now sigmagal
            PhEWP.machr = tab[idx]['Mach']
            PhEWP.msub = tab[idx]['LogMsub']
            PhEWP.mvir = tab[idx]['LogMvir']
            PhEWP.rvir = tab[idx]['Rvir']
            PhEWP.mcloud = tab[idx]['M_c']    

def draw(PhEWParticles, ax, fieldc, cmap='jet', \
         nskip=1, step=1, color_selection=True, \
         color_min=None, color_max=None, logcscale=False, alpha=0.2):
    if(color_min == None): color_min = min(fieldc)
    if(color_max == None): color_max = max(fieldc)
    phew.set_colors_phew_particles(PhEWParticles, fieldc, color_min, color_max, logscale=logcscale, cmap=plt.get_cmap(cmap))
    iskip = 0
    for PhEWP in PhEWParticles:
        if(iskip < nskip):
            iskip += 1
            continue
        flag = remove_spurious_particles(PhEWP.track)
        counter_spurious_particles[flag] += 1
        # flag = 1
        if(flag != 0): continue
        if(color_selection == True):
            if(fieldc == "Mvir"): value = PhEWP.mvir
            else: value = PhEWP.track[fieldc][0]
            if(not (color_min < value < color_max)):
                continue
        cs_a = sqrt(GAMMA * pc.k * PhEWP.track['T_a'] / (0.60 * pc.mh))
        a0 = PhEWP.track['atime'][0]
        # xarr = PhEWP.track['vrel'] / cs_a
        # yarr = PhEWP.track['T_c'] / PhEWP.track['T_a']
        xarr = PhEWP.track['Mach']
        yarr = PhEWP.track['M_cloud'] / MCLOUD
        # xarr = PhEWP.track['atime'] - a0
        # yarr = PhEWP.track['dr']

        # if(PhEWP.atimei > 0):
        #     ax.plot(PhEWP.machi, 1.0, "^", color=PhEWP.color, markersize=8)
        #     xarr[0], yarr[0] = PhEWP.machi, 1.0

        if(PhEWP.atimer > 0):
            ax.plot(PhEWP.machr, PhEWP.mcloud, "*", color=PhEWP.color, markersize=8)
        else:
            ax.plot(xarr[-1], yarr[-1], "+", color=PhEWP.color, markersize=6)            
        ax.plot(xarr[::step], yarr[::step], ".-", color=PhEWP.color, alpha=alpha, markersize=3)
        # ax.text(xarr[-1], yarr[-1], "%5.3f"%(PhEWP.track['atime'][-1]), color=PhEWP.color, fontsize=6)
        # ax.plot(xarr[0], yarr[0], "^", color=PhEWP.color, markersize=8)
        # ax.plot(xarr[-1], yarr[-1], "*", color=PhEWP.color, markersize=8)
        iskip = 0

Fields = phew.PhEWFields()
info_fields = Fields.get_field_info(["atime","Key","dr","dv","vrel","Mach","M_cloud","rho_a","Mstar","T_c","T_a"], verbose=True)
tab = phew.read_phew_fields(filename, info_fields)
# subtab = tab[tab['T_c'] < 1.e4]
pparts = phew.create_phew_particles(tab)

get_initwinds_and_rejoin_info(pparts, fphewsname)
# We want, Mvir, Mach (init), Mach (rejoin), f_cloud (rejoin)

fig = plt.figure(1, figsize=(6,7))
ax = fig.add_subplot(111)
# draw(pparts, ax, "T_a", nskip=10, logcscale=True, color_selection=True, color_min=2.e4, color_max=1.e7, alpha=0.3)

parts_to_show = select_particles(pparts, ntot=240)
draw(parts_to_show, ax, "Mvir", nskip=6, step=1, logcscale=False, color_selection=True, color_min=11.0, color_max=13.0, alpha=0.5, cmap=CMAP)
plt.text(0.6, 0.05, model, fontsize=12, transform = ax.transAxes)

fig.subplots_adjust(top=0.9, bottom=0.25)
axcbar = fig.add_axes([0.15,0.13,0.7,0.015])
norm1 = mpl.colors.Normalize(vmin=11.0, vmax=13.0)
cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.get_cmap(CMAP), norm=norm1, orientation="horizontal")
cdcbar.set_ticks([11, 11.5, 12, 12.5, 13])
cdcbar.set_ticklabels(["11","11.5","12","12.5","13"])
cdcbar.set_label("Mvir")

print (counter_spurious_particles)
ax.set_xlim(0.0, 20.0)
ax.set_ylim(0.0, 1.0)
ax.plot([1.0, 1.0], [0.0, 1.0], "k:") # vertical line
ax.plot([0.0, 20.0], [0.1, 0.1], "k:")
# ax.set_ylim(1.e-3, 2.0)
# ax.set_yscale("log")
ax.set_xlabel(r"$\mathcal{M}$")
ax.set_ylabel(r"$M_c/M_c(0)$")
plt.savefig('/home/shuiyao_umass_edu/figures/tmp.pdf')
plt.show()
