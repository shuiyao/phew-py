import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from pylab import setp
from astroconst import pc, ac
from scipy import interpolate

TMAX_THRESH = 5.5
PLOT_MODE = 0 # Cold/Hot/Wind
#PLOT_MODE = 1 # Metal
idxlist = [1,3]
#idxlist = [0,2]

def tcosmic(a):
    if(a == 0): return 0
    else: return 13.7 * a / sqrt(0.30/(a**3) + 0.70)

class acc_halo:
    def __init__(self, mh):
        self.mh = mh
        self.cold = array([0.]*100)
        self.hot = array([0.]*100)
        self.wind = array([0.]*100)
        self.windret = array([0.]*100)
        self.zwind = array([0.]*100)
        # self.zprimo = array([0.]*100)
        self.zcold = array([0.]*100)
        self.zhot = array([0.]*100)
        self.dtime = array([0.]*100)        

ztab, ttab = ioformat.rcol("/home/shuiyao/code/Python/tcosmic.dat", [1,2], linestart=1)
tcosmic = interpolate.interp1d(ztab, ttab)

ps = []        
simlist = [ \
        "p25n144rwlIII",\
        "p25n144phewIV"
        ]
fformat = [2, 2, 2]
#lines = [":","--", "-"]
lines = ["--","-", "-"]
lw = [1,1,1]
fig = plt.figure(1, figsize=(6,8))
zbins = linspace(0., 10., 100)
mhlist = ["105", "110", "115", "inf"]
axs = []
axs.append(fig.add_subplot(211))
axs.append(fig.add_subplot(212))
# Loop over each simulation:
for fidx in range(len(simlist)):
    acchalos = [acc_halo(10.5)]
    acchalos.append(acc_halo(11.0))
    acchalos.append(acc_halo(11.5))
    acchalos.append(acc_halo(12.0))
    for idx in range(4):
        fname = "/scratch/shuiyao/scidata/newwind/acc."+simlist[fidx]+"."+mhlist[idx]
        print "Reading: ", fname
        cold, hot, wind, windret, zcold, zhot, zwind = \
            ioformat.rcol(fname, [1,2,3,4,5,6,7])
        acchalos[idx].cold = cold
        acchalos[idx].zcold = zcold
        acchalos[idx].hot = hot
        acchalos[idx].zhot = zhot
        acchalos[idx].wind = wind
        acchalos[idx].zwind = zwind
        acchalos[idx].windret = windret

    axidx = 0
    for idx in idxlist:
        if(PLOT_MODE == 0): # Accretion Channels
            axs[axidx].plot(zbins, acchalos[idx].hot, color="red", linestyle=lines[fidx], linewidth=lw[fidx])
            p, = axs[axidx].plot(zbins, acchalos[idx].cold, color="blue", linestyle=lines[fidx], linewidth=lw[fidx])
            ps.append(p)
            axs[axidx].plot(zbins, acchalos[idx].wind, color="green", linestyle=lines[fidx], linewidth=lw[fidx])
            axs[axidx].plot(zbins, acchalos[idx].wind, color="green", linestyle=lines[fidx], linewidth=lw[fidx])
            axs[axidx].plot(zbins, acchalos[idx].windret, color="purple", linestyle=lines[fidx], linewidth=lw[fidx])
            axs[axidx].set_ylim(1.e-4, 1.0)        
        # plt.axis([0.,5.,1.e7,1.e11])
        if(PLOT_MODE == 1):
            axs[axidx].plot(zbins, array(acchalos[idx].zcold)/array(acchalos[idx].cold), color="blue", linestyle=lines[fidx], linewidth=lw[fidx])
            axs[axidx].plot(zbins, array(acchalos[idx].zhot)/array(acchalos[idx].hot), color="red", linestyle=lines[fidx], linewidth=lw[fidx])
            p, = axs[axidx].plot(zbins, array(acchalos[idx].zwind)/array(acchalos[idx].windret), color="green", linestyle=lines[fidx], linewidth=lw[fidx])
            ps.append(p)
            axs[axidx].set_ylim(5.e-6, 5.e-1)                
        axs[axidx].set_xlim(0., 5.)
        axs[axidx].set_yscale("log")
        axidx += 1

# plt.xscale("log")
plt.subplots_adjust(hspace=0)
setp(axs[0].get_xticklabels(),visible=False)
plt.xlabel("z")
#plt.ylabel("Accetion Rate")
# axs[0].text(0.1, 1.e10, "10.5 < Mfof < 11.0", fontsize=12)
# axs[1].text(0.1, 1.e10, "11.5 < Mfof", fontsize=12)
axs[0].text(0.1, 0.1, "10.5 < Mfof < 11.0", fontsize=12, transform=axs[0].transAxes)
axs[1].text(0.1, 0.1, "11.5 < Mfof", fontsize=12, transform=axs[1].transAxes)    
if(PLOT_MODE == 0):
    plt.ylabel("Acc Rate [Msolar/yr/Mpc^3]")
    axs[0].text(0.1, 0.26, "Cold (Primodial)", fontsize=10, transform=axs[0].transAxes, color="blue")
    axs[0].text(0.1, 0.22, "Hot (Primodial)", fontsize=10, transform=axs[0].transAxes, color="red")
    axs[0].text(0.1, 0.18, "Wind (All)", fontsize=10, transform=axs[0].transAxes, color="green")
    axs[0].text(0.1, 0.14, "Wind (Returned)", fontsize=10, transform=axs[0].transAxes, color="purple")
if(PLOT_MODE == 1):
    axs[0].text(0.1, 0.22, "wind", color="green", fontsize=12, transform=axs[0].transAxes)
    axs[0].text(0.1, 0.18, "hot", color="red", fontsize=12, transform=axs[0].transAxes)
    axs[0].text(0.1, 0.14, "cold", color="blue", fontsize=12, transform=axs[0].transAxes)
    plt.ylabel("Z_acc")

# plt.legend([ps[1],ps[2],ps[3]],["P25N144GW", "P25N144RWXI", "P25N144RWXX"], fontsize=12)
plt.legend([ps[1],ps[2]],simlist, fontsize=12)

# plt.title("0.45 < a < 0.5; 3/512 outputs")
plt.show()
        
