import matplotlib.pyplot as plt
import bin1d
from scipy import logspace, log10, array
import ioformat

import config_mpl
# import matplotlib as mpl
#import running_median
# mpl.rcParams['mathtext.default'] = "tt"
# mpl.rcParams['axes.labelsize'] = "large"

SCATTERPLOT = True

PLOT_ERB06 = True
PLOT_T04 = True

# 0: gadget3; 1: gizmo-g3wind-grackle; 6: l50n288
# 2: g3coolsf
# 1:  gizmo-g3wind-grackle blue
# 4:  phew purple
# 10  l25n144-phew-fkh100
# 11  l25n144-phew-zmixoff
# 12  l25n144-phew-condoff
# 13  l25n144-phew-m5kh100fs10
# 14  l25n144-phew-m4kh100fs10
# 15  l25n144-phew-m4kh50fs10
# 16  l25n144-phew-m4T3
# 17  l25n144-phew-m4kh100fs100
# 18  l25n144-phew-m5kh30fs10
# 19  l25n288-phew-m5kh30fs10
# SHOW_MODEL_LIST = [1, 13, 14, 15]
# SHOW_MODEL_LIST = [1, 12, 14, 15, 18, 19]
SHOW_MODEL_LIST = [1, 14, 18, 19]

#unit_m_l50 = 3469578.81574 * 1.e10
unit_m_l25 = 433697.351968 * 1.e10
#unit_m = 433697.351968 * 1.e10 * (12./25.) ** 3

from numpy import loadtxt
midx, models, clrs, lgds = loadtxt("models.dat", unpack=True, dtype='i8, S30, S20, S20')

texts = [ \
    "z=4.0 ", \
    "z=2.0", \
    "z=1.0", \
    "z=0.0" \
    ]
y_O = 0.009618
OHSOLAR = 8.69

ferb = open("../mzr/erb.2006.dat", "r")
Merb, Zerb, err = [], [], []
for line in ferb:
    Merb.append(float(line.split()[0]))
    Zerb.append(float(line.split()[2])-8.69)
    err.append(float(line.split()[3]))

def mzr_tremonti04(M):
    return -1.492 + 1.847*log10(M) - 0.08026*(log10(M)**2)

Zt04 = []
Mt04 = logspace(8.5, 11.5, 50)
for m in Mt04:
    Zt04.append(mzr_tremonti04(m)-8.69)

fig = plt.figure(1, figsize=(10,8))
plist, pflist = [], []
axs = []
for i in range(4):
    ax = plt.subplot(2,2,i+1)
    axs.append(ax)
    plt.vlines(log10(3.0e9), -1.5, 1.0,linestyles=":")
    plt.text(12.0, 0.6, texts[i], fontsize=16)
    ax.set_xlim(8.,13.)
    ax.set_ylim(-1.5,1.0)
YLABEL = "12+Log(O/H)"
axs[0].set_ylabel("[O/H]")
axs[2].set_ylabel("[O/H]")
axs[2].set_xlabel(r'$M_*/M_\odot$')
axs[3].set_xlabel(r'$M_*/M_\odot$')

def plot_mzr_model(mi, axs, shade=False, mass_unit=unit_m_l25):
    model = models[mi]
    clr = clrs[mi]
    mzrfile = [ \
                model+"/"+"mzr_033.txt", \
                model+"/"+"mzr_058.txt", \
                model+"/"+"mzr_078.txt", \
                model+"/"+"mzr_108.txt"]
    for i in range(len(mzrfile)):
        xline, yline = [], []
        f = open(mzrfile[i], "r")
        x, y = [], []
        if SCATTERPLOT == True:
            for line in f:
                spt = line.split()
                # ================================
                # The ACTUAL READING LINES
                # ================================
                x.append(float(spt[1])/0.7*mass_unit) # Gal Mass
                y.append(float(spt[3])) # idx, Mvir, C, <O>, Si, Fe
            for j in range(len(y)): # i is preserved!!
                if(y[j] > 0): y[j] = log10(y[j]/11.34) + 12. - 8.69

            x, y = bin1d.subsample(x, y, nonzero=True)
            s = bin1d.bin1d(x, y, nbins=40, logbins=True, bounds=0.95)
            xline, yline = s.cen, s.median
            # uline, lline = s.ubound, s.lbound
            xfline, uline, lline = [], [], []
            for idx in range(len(s.ubound)):
                if(-1.5 < s.ubound[idx] < 1.0): 
                    if(-1.5 < s.lbound[idx] < 1.0):
                        lline.append(s.lbound[idx])
                        uline.append(s.ubound[idx])
                        xfline.append(log10(s.cen[idx]))
        f.close()
        xline = log10(xline)
        p, = axs[i].plot(xline, yline, "-", color=clr) # Draw line for all models
        plist.append(p)
        if(shade == True):
            pf, = plt.fill(xfline+xfline[::-1], uline+lline[::-1], color=colors_scatter[i], alpha=0.2)
            pflist.append(pf)
    

def convert_pp04n2_t04(Z): # Z = 12 + log(O/H)
    # Kewley & Ellison 2008. Table 3
    a, b, c, d = -1661.9380, 585.17650, -68.471750, 2.6766690
    # if(8.05 < Z < 8.9):
    return a + b * Z + c * Z * Z + d * Z ** 3

for mi in SHOW_MODEL_LIST:
    mass_unit = unit_m_l25
    if(mi == 19): mass_unit = unit_m_l25
    plot_mzr_model(mi, axs, mass_unit=mass_unit)

if PLOT_ERB06 == True:    
    ferb = open("../mzr/erb.2006.dat", "r")
    Merb, Zerb, err = [], [], []
    Merr = []
    for line in ferb:
        Mcen = float(line.split()[0])
        Merb.append(log10(Mcen))
        Merr.append(log10(Mcen+float(line.split()[1]))-log10(Mcen))
        Zerb.append(float(line.split()[2])-0.0)
        err.append(float(line.split()[3]))
    # ax.plot(Merb, Zerb, "ko")
    Zmid = convert_pp04n2_t04(array(Zerb))
    Zupper = convert_pp04n2_t04(array(Zerb)+array(err))
    Zlower = convert_pp04n2_t04(array(Zerb)-array(err))
    err1 = Zupper - Zmid
    err2 = Zmid - Zlower
    axs[1].errorbar(Merb, Zmid-OHSOLAR, xerr=Merr, yerr=[err1, err2], color="black", fmt='o')
    axs[1].errorbar(Merb[0], Zmid[0]-OHSOLAR, xerr=Merr[0], yerr=0.05, uplims=True, color="black", fmt='o')
    # axs[1].plot(Merb, Zmid-OHSOLAR-0.3, color="purple", linestyle="-", linewidth=2)
    axs[1].text(0.02, 0.02, "Erb+ (2006)", transform=axs[1].transAxes, fontsize=12, color="black")
if PLOT_T04 ==True:
    ax = axs[3]
    Mt04, Z16t04, Z50t04, Z84t04 = ioformat.rcol("../mzr/tremonti04.dat", [0,2,3,4], linestart=1)
    Z16t04 = array(Z16t04) - OHSOLAR
    Z50t04 = array(Z50t04) - OHSOLAR
    Z84t04 = array(Z84t04) - OHSOLAR
    p, = ax.plot(Mt04, Z50t04, color="black", linestyle="-")
    # ax.plot(Mt04, Z50t04-0.3, color="purple", linestyle="-", linewidth=2)
    ax.plot(Mt04, Z16t04, color="black", linestyle="--")
    ax.plot(Mt04, Z84t04, color="black", linestyle="--")
    plist.append(p)
    axs[3].text(0.02, 0.02, "Tremonti+ (2004)", transform=axs[3].transAxes, fontsize=12, color="black")

from pltastro import legend
lgd = legend.legend(axs[0])
lgd.loc="lower right"
for mi in SHOW_MODEL_LIST:
    lgd.addLine((lgds[mi], clrs[mi], "-", 1))
lgd.draw()
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
