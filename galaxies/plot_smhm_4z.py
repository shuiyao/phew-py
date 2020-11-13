import ioformat
import matplotlib.pyplot as plt
import bin1d
from scipy import logspace, log10, append, linspace, array, loadtxt
import matplotlib as mpl
import behroozi
from pylab import setp
from matplotlib import gridspec

mpl.rcParams['mathtext.default'] = "tt"
mpl.rcParams['axes.labelsize'] = "large"

PLOT_MASS_RATIO = True

# SHOW_MODEL_LIST = [1, 3, 13, 14]
# SHOW_MODEL_LIST = [1, 3, 10, 11, 13, 14] # sw

# SPH method. PhEWOff x PhEW
# SHOW_MODEL_LIST = [3, 7, 6, 18]

# PHEWs
# SHOW_MODEL_LIST = [6, 19]
#SHOW_MODEL_LIST = [1, 2, 3, 4, 5, 6]
SHOW_MODEL_LIST = [6, 13, 14]

zstr = ["036", "005", "003", "001"]
zstr2 = ["108", "078", "058", "033"]
zs = [0.0, 1.0, 2.0, 4.0]

fbase = "/scratch/shuiyao/sci/PHEW_TEST/"
FILLED = False
#HALOMASS_CORRECTION = True
#modelname = "l25n144-phewoff"
#modelname = "l25n144-gizmo-mufasa"
#modelname = "l25n144-mufasa"

HALOMASS_CORRECTION = True
# modelname = "l25n144-phewoff-nozload"
# modelname = "l25n144-gizmo-mufasa"
modelname = "l25n144-mufasa"

midx, models, clrs, lgds = loadtxt("models.dat", unpack=True, dtype='i8, S30, S20, S30')

flist = [\
          fbase + modelname + "/" + "smhm_" + zstr2[0]+".txt",
          fbase + modelname + "/" + "smhm_" + zstr2[1]+".txt",
          fbase + modelname + "/" + "smhm_" + zstr2[2]+".txt",
          fbase + modelname + "/" + "smhm_" + zstr2[3]+".txt"          
         ]

MLIM = 2.*5.8e9 # 128 Msph
MLIM_PLOT = 2.9e9 # 32 Msph
MLIM = 2.*5.8e9 / 8.# 128 Msph
MLIM_PLOT = 2.9e9 / 8.# 32 Msph
PLOT_MOSTER = 0

colors = [ \
    "black",
    "red", \
    "blue", \
    "purple", \
    "black", \
    "red", \
    "blue", \
    "purple", \                                                       
]

fig = plt.figure(1, figsize=(8,8))
axs = []
for i in range(4):
    axs.append(fig.add_subplot(2,2,i+1))

def plot_smhm_model(mi, axs, shade=False, mhcorr=1.0):
    modelname = models[mi]
    clr = clrs[mi]
    smhmlist = [\
             fbase + modelname + "/" + "smhm_" + zstr2[0]+".txt",
             fbase + modelname + "/" + "smhm_" + zstr2[1]+".txt",
             fbase + modelname + "/" + "smhm_" + zstr2[2]+".txt",
             fbase + modelname + "/" + "smhm_" + zstr2[3]+".txt"          
    ]
    for i in range(len(smhmlist)):
        # ms, mvir, flag, mh = ioformat.rcol(smhmlist[i], [0,1,2,3], linestart=1)
        ms, mh, flag = ioformat.rcol(smhmlist[i], [0,1,2], linestart=1)
        x, y = [], []
        mh = array(mh) / mhcorr
        for j in range(len(flag)):
            if(flag[j] == 1.0): # Only plot the central
                x.append(mh[j])
                if(PLOT_MASS_RATIO == True): y.append(ms[j]/mh[j])
                else: y.append(ms[j])
        # plt.plot(x, y, "b.", alpha=0.2)
        x, y = bin1d.subsample(x, y, nonzero=True)
        s = bin1d.bin1d(x, y, nbins=40, logbins=True, bounds=0.95)
        xline, yline = log10(s.cen), log10(s.median)
        uline, lline = log10(s.ubound), log10(s.lbound)
        axs[i].plot(xline, yline, "-", color=clr, linewidth=2)
        maxidx = -4
        if(shade == True):
            axs[i].fill(append(xline[:maxidx],xline[maxidx-1::maxidx+1]), append(uline[:maxidx],lline[maxidx-1::maxidx+1]), color=colors[i], alpha=0.2)

for mi in SHOW_MODEL_LIST:
    plot_smhm_model(mi, axs, shade=False, mhcorr=1.0) 
# plot_smhm_model(0, axs, shade=False, mhcorr=8.0) # l25n144-gadget3
# plot_smhm_model(1, axs, shade=False, mhcorr=1.0) # l25n144-gizmo-mufasa
# plot_smhm_model(5, axs, shade=False, mhcorr=1.0) # l25n144-nosniifac
# plot_smhm_model(12, axs, shade=False, mhcorr=1.0/8.0) # l25n144-mufasa-dust
# plot_smhm_model(13, axs, shade=False, mhcorr=1.0/8.0) # l25n144-dust

txt = ["0", "1", "2", "4"]    
for i in range(4):
    plt.hlines(log10(MLIM), 11., 13.5, color="k", linestyles="--")
    axs[i].text(10.7, -1.2, "z = "+txt[i])    
    if(PLOT_MASS_RATIO == True):
        axs[i].set_xlim(10.5, 13.5)
        axs[i].set_ylim(-4, -1)            
    else:
        plt.axis([10.5, 13.5, log10(MLIM_PLOT), 12.0])

if(PLOT_MASS_RATIO == True):
    for i in range(4):
        mh, dm, e1, e2 = behroozi.read_smmr(zs[i])
        axs[i].errorbar(mh, dm, yerr=[e1, e2], color="black", fmt='o')        
        axs[i].plot(mh, dm, "-", color="black")
        # if(i == 1):
        #     mh, ms = ioformat.rcol("../REFERENCES/moster_2013_2018/moster18_z1.dat", [0, 1])
        #     dm = array(ms) - array(mh)
        #     axs[i].plot(mh, dm, "-", color="lightgrey")

#plt.legend(["PESPH-HM12-AC", "P50N576GW"])
# plt.legend([p1, p3, p4, p2, p5], ["DESPH", "PESPH-HM12-AC", "PESPH-GW-5xEsn", "DESPH", "PESPH-GW-5xEsn"], fontsize=12, loc=2)
# plt.legend([p1, p3, p4, p2, p5], ["DESPH", "PESPH", "PESPH-GW", "DESPH", "PESPH-GW"], fontsize=20, loc=2)
if(PLOT_MASS_RATIO == False):
    axs[0].set_ylabel(r'$Log(M*/M_\odot)$')
    axs[2].set_ylabel(r'$Log(M*/M_\odot)$')
else:
    axs[0].set_ylabel(r'$Log(M*/M_h)$')
    axs[2].set_ylabel(r'$Log(M*/M_h)$')
axs[2].set_xlabel("$Log(M_h/M_\odot)$")
axs[3].set_xlabel("$Log(M_h/M_\odot)$")
fig.subplots_adjust(hspace=0.0, wspace=0.0)
setp(axs[0].get_xticklabels(),visible=False)
setp(axs[1].get_xticklabels(),visible=False)
setp(axs[1].get_yticklabels(),visible=False)
setp(axs[3].get_yticklabels(),visible=False)
axs[0].yaxis.set_ticks(linspace(-4.0, -1.0, 4))
axs[2].yaxis.set_ticks(linspace(-4.0, -2.0, 3))

from pltastro import legend
lgd = legend.legend(axs[0])
lgd.loc="lower right"
for mi in SHOW_MODEL_LIST:
    lgd.addLine((lgds[mi], clrs[mi], "-", 1))
lgd.draw()
# plt.savefig("smms.pdf")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
