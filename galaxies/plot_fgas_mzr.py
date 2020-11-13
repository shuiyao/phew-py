import matplotlib.pyplot as plt
import bin1d
import matplotlib as mpl
from scipy import log10, array, concatenate, linspace, isfinite
import ioformat
#import running_median
# mpl.rcParams['mathtext.default'] = "tt"
# mpl.rcParams['axes.labelsize'] = "large"
import config_mpl

PLOT_T04 = True
PLOT_ERB06 = False
PLOT_SANDERS15 = True
PLOT_ZAHID13 = False
PLOT_P14 = True

STELLARMASS = True
HALOMASS = False
#OHSOLAR = 8.69
OHSOLAR = 0.0
unit_m1 = 3469578.81574 * 1.e10 / 8.
unit_m2 = 3469578.81574 * 1.e10 / 8.
axs = []

# lgds = ["ezwDESPH", "RefHres", "Ref", "RefHres"]
lgds = ["M5", "M4", "Gadget3", "M4"]
files = []
# models = ["l50n288-phew-m5", "l50n288-phew-m4"]
models = ["l25n288-phew-m5", "l25n288-phew-m4"]
zstrs = ["058", "108"]
for zstr in zstrs:
    for model in models:
        fname = model + "/fgas_mstar_" + zstr + ".txt"
        print fname
        files.append(fname)
    for model in models:
        fname = model + "/mzr_" + zstr + ".txt"
        print fname
        files.append(fname)

p50n288zwfiles = [\
    "../fgas/fgas_mstar_p50n288zw_058.txt", \          
    "../mzr/massz_p50n288zw_058.txt", \          
    "../fgas/fgas_mstar_p50n288zw_108.txt", \          
    "../mzr/massz_p50n288zw_108.txt", \                    
]
colors = [ \
    "blue", \
    "black", \
    "blue", \
    "black", \
    "blue", \
    "black", \
    "blue", \
    "black" \
    ]
colors_scatter = [ \
    "blue", \
    "black", \
    "blue", \
    "black", \
    "blue", \
    "black", \
    "blue", \
    "black", \
    ]
linestyles = [ \
    "-", \
    "-", \
    "-", \
    "-", \
    "-", \
    "-", \
    "-", \
    "-" \
    ]
texts = [ \
    "z=2.0 ", \
    "z=2.0", \
    "z=0.0", \
    "z=0.0" \
    ]

def frac_red(mstar, fgas):
    ntot, nred, nredmassive = 0, 0, 0
    c = zip(mstar, fgas)
    c.sort()
    mstar, fgas = zip(*c)
    ntot = len(mstar)
    for i in range(ntot):
        if(fgas[i] == 0): nred += 1
    print "red fraction: ", float(nred)/float(ntot)
    nred = 0.0
    for i in range(int(ntot*0.01)):
        if(fgas[-i-1] == 0): nred += 1
    print "red fraction (massive 1%): ", float(nred)/(float(ntot) * 0.01)
        
fig = plt.figure(1, figsize=(9,8))
fig.set_rasterized(True)
plist, pflist = [], []
for i in range(len(files)):
    panel = int(i/2)+1
    if panel == i/2.+1.: # i = 0, 2, 4, 6
        unit_m = unit_m1
        ax = plt.subplot(2,2,panel)
        axs.append(ax)
        xline, yline = [], []
    else:
        unit_m = unit_m2
    f = open(files[i], "r")
    x, y = [], []
    xsat, ysat = [], []

    if(panel == 2 or panel == 4): # MZR
        print "i=%d, panel=%d, MZR" % (i, panel)
    else:
        print "i=%d, panel=%d, Fgas" % (i, panel)            
    
    for line in f:
        spt = line.split()
        # ================================
        # The ACTUAL READING LINES
        # ================================
        y.append(float(spt[3]))
        if(panel == 2 or panel == 4): # MZR
            x.append(float(spt[1])/0.7*unit_m) # Stellar Mass
        else:
            # x.append(log10(float(spt[2])/0.7*unit_m)) # Stellar Mass
            x.append(float(spt[2])/0.7*unit_m) # Stellar Mass            
    if(panel == 2 or panel == 4): # MZR
        for j in range(len(y)): # i is preserved!!
            if(y[j] > 0): y[j] = log10(y[j]/11.34) + 12. - OHSOLAR
    else:
        frac_red(x, y)
    x, y = bin1d.subsample(x, y, nonzero=True)
    s = bin1d.bin1d(x, y, nbins=40, logbins=True, bounds=0.84)
    xline, yline = s.cen, s.median
    # xline, yline = s.value, s.median
    uline, lline = s.ubound, s.lbound
    f.close()
    p, = ax.plot(log10(array(xline)), yline, linestyles[i], color=colors[i]) # Draw line for all models
    plist.append(p)
    # if(panel == 3): axs[0].plot(log10(array(xline)), yline, "--", color=colors[i])
    lcount = i-int(i/2)*2
    # plt.vlines(3.e9, 0., 1.,linestyles=":")
    plt.text(0.05, 0.90, texts[panel-1], transform=ax.transAxes, fontsize=12)
    ax.set_xlim(9.,12.5)
    ax.set_xticks([9.0, 10.0, 11.0, 12.0])
    if panel == 1 or panel == 3:
        ax.set_ylim(0.,1.)
        ax.set_ylabel(r'$f_{gas}$')
    if panel == 2 or panel == 4:
        ax.set_ylim(8.0,9.5)
        ax.set_ylabel(r'12+log([O/H])')
    if panel == 4 or panel == 3:
        ax.set_xlabel(r'$Log(M_*/M_\odot)$')
    xfill, ufill, lfill = [], [], []
    xline = log10(array(xline))
    for j in range(len(xline)):
        if(isfinite(uline[j])):
            if(isfinite(lline[j])):
                xfill.append(xline[j])
                ufill.append(uline[j])
                lfill.append(lline[j])                
    # pf, = plt.fill(log10(array(xline+xline[::-1])), uline+lline[::-1], color=colors_scatter[i], alpha=0.4)
    if(i in [1, 3, 5, 7]):
        pf, = plt.fill(xfill+xfill[::-1], ufill+lfill[::-1], color=colors_scatter[i], alpha=0.3)
        pflist.append(pf)

# Plot the dashed lines        
for i in range(len(p50n288zwfiles)):
    f = open(p50n288zwfiles[i], "r")
    x, y = [], []
    xline, yline = [], []    
    for line in f:
        spt = line.split()
        y.append(float(spt[3]))
        if(i == 1 or i == 3): # Fgas
            x.append(float(spt[1])/0.7*unit_m) # Stellar Mass
        else:
            x.append(float(spt[2])/0.7*unit_m) # Stellar Mass            
    if(i == 1 or i == 3): # MZR
        for j in range(len(y)): # i is preserved!!
            if(y[j] > 0): y[j] = log10(y[j]/11.34) + 12. - OHSOLAR
    else:
        frac_red(x, y)
    x, y = bin1d.subsample(x, y, nonzero=True)
    s = bin1d.bin1d(x, y, nbins=40, logbins=True, bounds=0.84)
    xline, yline = s.cen, s.median
    uline, lline = s.ubound, s.lbound
    f.close()
    p, = axs[i].plot(log10(array(xline)), yline, "-", color="magenta") # Draw line for all models
    if(i == 0):
        plist.append(p)

# axs[2].legend([plist[0], plist[1], plist[-1], pflist[0], pflist[1]], ["ezwDESPH","RefHres","Ref", "ezwDESPH","RefHres"], prop={'size':12}, loc=1)
# axs[2].legend([plist[0], plist[1], plist[-1], pflist[0]], ["ezwDESPH","RefHres","Ref", "RefHres"], prop={'size':12}, loc=1)
axs[2].legend([plist[0], plist[1], plist[-1], pflist[0]], lgds, prop={'size':12}, loc=1)

if PLOT_T04 ==True:
    Mt04, Z16t04, Z50t04, Z84t04 = ioformat.rcol("../mzr/tremonti04.dat", [0,2,3,4], linestart=1)
    Z16t04 = array(Z16t04) - OHSOLAR
    Z50t04 = array(Z50t04) - OHSOLAR
    Z84t04 = array(Z84t04) - OHSOLAR
    p, = ax.plot(Mt04, Z50t04, color="green", linestyle="-")
    ax.plot(Mt04, Z16t04, color="green", linestyle="--")
    ax.plot(Mt04, Z84t04, color="green", linestyle="--")
    plist.append(p)
    axs[3].text(0.02, 0.02, "Tremonti+ (2004)", transform=axs[3].transAxes, fontsize=12, color="green")    

if PLOT_P14 == True:
    # Mp14, F50p14, F16p14, F84p14 = ioformat.rcol("../fgas/peeples14.dat", [0,1,2,3], linestart=5)
    from numpy import loadtxt
    Mp14, F50p14, F16p14, F84p14 = loadtxt("../fgas/peeples14.dat", unpack=True)    
    F50p14 = array(F50p14) / (1. + array(F50p14))
    F16p14 = array(F16p14) / (1. + array(F16p14))
    F84p14 = array(F84p14) / (1. + array(F84p14))
    err1 = F16p14 - F50p14
    err2 = F50p14 - F84p14    
    axs[2].errorbar(Mp14, F50p14, [err1, err2], color="green", fmt='o')
    axs[2].text(0.02, 0.02, "Peeples+ (2014)", transform=axs[2].transAxes, fontsize=12, color="green")    

def convert_pp04n2_t04(Z): # Z = 12 + log(O/H)
    # Kewley & Ellison 2008. Table 3
    a, b, c, d = -1661.9380, 585.17650, -68.471750, 2.6766690
    # if(8.05 < Z < 8.9):
    return a + b * Z + c * Z * Z + d * Z ** 3

if PLOT_ERB06 == True:    
    ferb = open("../mzr/erb.2006.dat", "r")
    Merb, Zerb, err = [], [], []
    Merr = []
    for line in ferb:
        Mcen = float(line.split()[0])
        Merb.append(log10(Mcen))
        Merr.append(log10(Mcen+float(line.split()[1]))-log10(Mcen))
        Zerb.append(float(line.split()[2])-OHSOLAR)
        err.append(float(line.split()[3]))
    # ax.plot(Merb, Zerb, "ko")
    Zmid = convert_pp04n2_t04(array(Zerb))
    Zupper = convert_pp04n2_t04(array(Zerb)+array(err))
    Zlower = convert_pp04n2_t04(array(Zerb)-array(err))
    err1 = Zupper - Zmid
    err2 = Zmid - Zlower
    axs[1].errorbar(Merb, Zmid, xerr=Merr, yerr=[err1, err2], color="green", fmt='o')
    axs[1].errorbar(Merb[0], Zmid[0], xerr=Merr[0], yerr=0.05, uplims=True, color="green", fmt='o')
    axs[1].text(0.02, 0.02, "Erb+ (2006)", transform=axs[1].transAxes, fontsize=12, color="green")

if PLOT_SANDERS15 == True:    
    fsanders = "../mzr/sanders15n2.dat"
    M0, M1, M2, Z0, Z1, Z2 = \
        ioformat.rcol(fsanders, [0,1,2,3,4,5], linestart=3)    
    Zmid = convert_pp04n2_t04(array(Z0))
    Zupper = convert_pp04n2_t04(array(Z0)+array(Z2))
    Zlower = convert_pp04n2_t04(array(Z0)+array(Z1))
    Merr1 = array(M0) - array(M1)
    Merr2 = array(M2) - array(M0)
    Zerr1 = Zmid - Zlower
    Zerr2 = Zupper - Zmid    
    axs[1].errorbar(M0, Zmid, xerr=[Merr1, Merr2], yerr=[Zerr1, Zerr2], color="green", fmt='o')
    axs[1].text(0.02, 0.02, "Sanders+ (2015)", transform=axs[1].transAxes, fontsize=12, color="green")

# Zahid et al., 2013, ApJL, 771, L19    
def zahid13(M, kind):
    if(kind == "SDSS"):
        Z0, M0, g = 9.121, 8.999, 0.85
    if(kind == "E06"):
        Z0, M0, g = 9.06, 9.7, 0.6
    return Z0 - log10(1. + (10.**(M - M0))**(-g))

if PLOT_ZAHID13 == True:
    Mx = linspace(9.0, 10.5, 40)
    Zy = zahid13(Mx, kind="E06")
    axs[1].plot(Mx, Zy, "r--")
    Zy = zahid13(Mx, kind="SDSS")
    axs[3].plot(Mx, Zy, "-", color="darkblue")
    
fig.subplots_adjust(wspace=0.40, hspace=0.15)
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
