from myinit import *
import matplotlib.pyplot as plt
from scipy import linspace, exp, log10, log, array, logspace
from scipy.integrate import quad
from scipy.special import gamma
from matplotlib import gridspec
from pylab import setp
import ioformat

# OrigGasMass = 0.00653561 - 9.3e7 Msolar
# 64 Particles = 4.1827904e9M_solar

YMIN, YMAX = -5.0, -1.0
#YMIN, YMAX = -5.4, 0.6

SHOW_MODEL_LIST = [1, 2, 3, 4, 5, 6, 8]
lw = [1, 2, 1, 1, 2, 2, 2]
lstyles = ["--", "--", "-", "-", "-", "-", "-"]
#SHOW_MODEL_LIST = [0, 1, 4, 8]
fbase = DIRS['SCIDATA']

#P50N288
MLIM_LRES = 2.*5.8e9 # 128 Msph
MLIM_PLOT = 2.9e9 / 2. # 32 Msph
# P25N288
MLIM_HRES = 1.49e9 # 128 Msph
MLIM_PLOT = 1.49e9 / 4. # 32 Msph

MLIM_LRES = log10(MLIM_LRES)
MLIM_HRES = log10(MLIM_HRES)
MLIM_PLOT = log10(MLIM_PLOT)

PLOT_MUFASA = False
PLOT_EAGLE = False
SAVE_FIGURE = False
PLOT_BERNARDI13 = False

from numpy import loadtxt
midx, models, clrs, lgds = loadtxt("models.dat", unpack=True, dtype='i8, U30, U20, U20')

titlestr = "P50N288, Baldry+12, Tomczak+14, EAGLE"
fig = plt.figure(1, figsize=(8,8))
axs = []
for i in range(4):
    axs.append(fig.add_subplot(2,2,i+1))

lines_z = []
lines_data = []
lines_runs = []
zstr2 = ["108", "078", "058", "033"]

# colors = ["black","black","black","black"]
# colors = ["black", "red", "blue", "purple"]
colors = ["black", "darkgray", "grey", "lightgrey"]

def read_gsmfs(fname):
    f = open(fname, "r")
    f.readline()
    Nbins = int(f.readline().split()[0])
    x, y = [], []
    for j in range(Nbins):
        line = f.readline()
        phi = float(line.split()[1])*0.7**3
        if(phi != 0):
            x.append(log10(float(line.split()[0])/0.7))
            y.append(log10(phi))
    return x, y

# -------------------------------- DATA --------------------------------
def plot_gsmf_model(mi, axs, noplotlist=[], lw=1, lstyle="-"):
    modelname = models[mi]
    clr = clrs[mi]
    gsmfs = [ \
              fbase + modelname + "/" + "gsmf_" + zstr2[0]+".txt",
              fbase + modelname + "/" + "gsmf_" + zstr2[1]+".txt",
              fbase + modelname + "/" + "gsmf_" + zstr2[2]+".txt",
              fbase + modelname + "/" + "gsmf_" + zstr2[3]+".txt"]
    for i in range(len(gsmfs)):
        if(i in noplotlist): continue
        print (gsmfs[i])
        ms, phi = read_gsmfs(gsmfs[i])
        line, = axs[i].plot(ms, phi, linestyle=lstyle, color=clr, linewidth=lw)

for i, mi in enumerate(SHOW_MODEL_LIST):
    plot_gsmf_model(mi, axs, lw=lw[i], lstyle=lstyles[i])

def schechter(M, pars):
    m_internal = M / 10.**(pars[0])
    c1 = exp(-m_internal)
    c2 = pars[1] * m_internal ** (pars[2])
    c3 = pars[3] * m_internal ** (pars[4])
    return c1 * (c2 + c3) / 1.e3# / 0.7**3

pars_baldry12 = [10.648, 4.26, -0.46, 0.58, -1.58]
pars_ilbert11_z0 = [10.88, 1.68, -0.69, 0.77, -1.42]
pars_ilbert11_z1 = [10.87, 2.03, -0.52, 0.29, -1.62]
pars_ilbert11_z2 = [10.74, 0.88, -0.24, 0.33, -1.6]
pars_t14_z0 = [10.78, 2.88, -0.98, 0.05, -1.9]
pars_t14_z1 = [10.54, 1.90, 0.30, 0.68, -1.45]
pars_t14_z2 = [10.69, 0.15, 1.03, 0.55, -1.33]

gsmf_baldry12 = "../../REFERENCES/Baldry12/gsmf_baldry12_z0.dat"
ms, phi, err = ioformat.rcol(gsmf_baldry12, [0,2,3], linestart=3)
phi = array(phi) / 1.e3
err= array(err) / 1.e3
lower = log10(phi) - log10(phi - err)
upper = log10(phi + err) - log10(phi)
phi = log10(phi)
line = axs[0].errorbar(ms, phi, yerr=[lower, upper], color="black", fmt='o')
lines_data.append(line) # Baldry +12

# Based on Bernardi.2013.MNRAS.436.697
# cmodel
def bernardi13(M, fitmodel="sersic"):
    if(fitmodel == "cmodel"):
        phi1, phi2 = 0.766e-2, 0.557e-2
        M1, M2 = 4.103e8, 4.7802e9
        a, b, c = 1.764, 0.384, 0.053
    if(fitmodel == "sersic"):    
        phi1, phi2 = 1.040e-2, 0.675e-2
        M1, M2 = 0.0094e9, 2.7031e9
        a, b, c = 1.665, 0.255, 0.296
    F1, F2 = M / M1, M / M2
    y = phi1 * b / gamma(a/b) * (F1) ** a * exp(-F1**b)
    y += phi2 * (F2) ** c * exp(-F2)
    return log10(log(10.)*y)

if(PLOT_BERNARDI13 == True):
    ms, phi = logspace(9., 13., 50), []
    for mi in ms: phi.append(bernardi13(mi, fitmodel="cmodel"))
    line_b13a, = axs[0].plot(log10(ms), phi, ":", color="green", linewidth=2)
    ms, phi = logspace(9., 13., 50), []
    for mi in ms: phi.append(bernardi13(mi, fitmodel="sersic"))
    line_b13b, = axs[0].plot(log10(ms), phi, "--", color="green", linewidth=2)

gsmfs_tomczak14 = [ \
    "../../REFERENCES/tomczak14/t14_z0.2_z0.5.dat",\
    "../../REFERENCES/tomczak14/t14c_z0.75_z1.25.dat",\
    "../../REFERENCES/tomczak14/t14c_z1.5_z2.5.dat",\
]
# colors = ["green", "red", "blue"]
axlst = [0,1,2]
for i in range(len(gsmfs_tomczak14))[:3]:
    ms, phi, e1, e2 = ioformat.rcol(gsmfs_tomczak14[i], [0,1,2,3], linestart=3)
    # phi = array(phi) + log10(1.2 * 1.13)
    # lower = 10.**(array(phi)) - 10.**(array(phi) + array(e2))
    # upper = 10.**(array(phi) + array(e1)) - 10.**(array(phi))
    # phi0 = 10.**(array(phi))
    lower = -array(e2)
    upper = array(e1)
    if(i > 0):
        line = axs[axlst[i]].errorbar(ms, phi, yerr=[lower, upper], color=colors[i], fmt='o')
    # plt.plot(ms, phi0, "o", color=colors[i])
        lines_data.append(line) # Tomczak +14

gsmf_song16 = "../../REFERENCES/song16/gsmf_song16_z4.dat"
ms, phi, e1, e2 = ioformat.rcol(gsmf_song16, [0,1,2,3], linestart=3)
lower = -array(e2)
upper = array(e1)
line = axs[3].errorbar(ms, phi, yerr=[lower, upper], color="grey", fmt='^')
lines_data.append(line) # Song +16

# MUFASA
def read_mufasa(fname):
    ms, phi = [], []
    f = open(fname, "r")
    f.readline()
    f.readline()
    ndata = int(f.readline().split()[0])
    for i in range(ndata): ms.append(float(f.readline().split()[0]))
    for i in range(ndata): phi.append(float(f.readline().split()[0]))
    f.close()
    return ms, phi

gsmfs_mufasa = [\
    "../../REFERENCES/MUFASA/gsmf_mufasa_z0.dat",\
    "../../REFERENCES/MUFASA/gsmf_mufasa_z1.dat",\
    "../../REFERENCES/MUFASA/gsmf_mufasa_z2.dat",\
    "../../REFERENCES/MUFASA/gsmf_mufasa_z4.dat"\                
]
# colors = ["black", "red", "blue", "purple"]
if(PLOT_MUFASA == True):
    for i in range(len(gsmfs_mufasa)):
        ms, phi = read_mufasa(gsmfs_mufasa[i])
        line, = axs[i].plot(ms, phi, "-", color=colors[i], linewidth=1)
        if(i == 0): lines_runs.append(line) # MUFASA
        
if(PLOT_EAGLE == True):
    gsmfs_eagle = [ \
        "../../REFERENCES/EAGLE/gsmf/Ref_gsmf_z0p1.txt",\
        "../../REFERENCES/EAGLE/gsmf/Ref_gsmf_z1p0.txt",\
        "../../REFERENCES/EAGLE/gsmf/Ref_gsmf_z2p0.txt",\
        "../../REFERENCES/EAGLE/gsmf/Ref_gsmf_z4p0.txt"
    ]
    # colors = ["black", "red", "blue", "purple"]
    for i in range(len(gsmfs_eagle)):
        ms, phi = ioformat.rcol(gsmfs_eagle[i], [0,2], linestart=2)
        ms = log10(array(ms))
        phi = log10(array(phi))
        line, = axs[i].plot(ms, phi, ":", color=colors[i])
        if(i == 0): lines_runs.append(line) # EAGLE

txt = ["0", "1", "2", "4"]    
for i in range(4):
    axs[i].set_xlim(MLIM_PLOT, 13.0)
    axs[i].set_ylim(YMIN, YMAX)
    axs[i].text(12.0, -1.5, "z = "+txt[i])
    axs[i].plot([MLIM_LRES, MLIM_LRES], [YMIN, YMAX], "k:")
    axs[i].plot([MLIM_HRES, MLIM_HRES], [YMIN, YMAX], "k--")    
#l1 = plt.legend(lines, ["z = 0", "z = 1", "z = 2"], loc=3, fontsize=12)
# l1 = axs[0].legend(lines_runs, [runname,"EAGLE","P50-N576"], loc=3, fontsize=12)
# l1 = axs[0].legend([line_b13a, line_b13b], ["Dave +16", "EAGLE", "B13, cmodel", "B13, sersic"], loc=3, fontsize=12)
l2 = axs[1].legend(lines_data, ["Baldry +12", "Tomczak +14", "Tomczak +14", "Song +16"], loc=3, fontsize=12)
# plt.gca().add_artist(l1)
axs[0].set_ylabel(r'$\Phi [Mpc^{-3}/Log(M)]$', fontsize=12)
axs[2].set_ylabel(r'$\Phi [Mpc^{-3}/Log(M)]$', fontsize=12)
axs[2].set_xlabel("Log(M*)")
axs[3].set_xlabel("Log(M*)")
fig.subplots_adjust(hspace=0.0, wspace=0.0)
setp(axs[0].get_xticklabels(),visible=False)
setp(axs[1].get_xticklabels(),visible=False)
setp(axs[1].get_yticklabels(),visible=False)
setp(axs[3].get_yticklabels(),visible=False)
axs[0].yaxis.set_ticks(linspace(-5.0, -1.0, 5))
axs[2].yaxis.set_ticks(linspace(-5.0, -2.0, 4))

from pltastro import legend
lgd = legend.legend(axs[3])
lgd.loc="lower right"
lgd.fontsize = 8
for i, mi in enumerate(SHOW_MODEL_LIST):
    lgd.addLine((lgds[mi], clrs[mi], lstyles[i], lw[i]))
lgd.draw()

# plt.title(titlestr)
plt.savefig(DIRS['FIGURE']+"tmp.pdf")

plt.show()

