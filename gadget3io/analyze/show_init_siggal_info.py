from myinit import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic
import matplotlib.pyplot as plt
import pandas as pd

# Pre-requisite: $SCI/phew-py/select/select_winds_dzfix.sh $modelname

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 40
SIGMAX = 150.0
sigbins = linspace(0., SIGMAX, nbins+1)
sigcen = 0.5 * (sigbins[1:] + sigbins[:-1])
dsig = SIGMAX / nbins

ZSOLAR = log10(0.0122)

def find_siggal_bin(sig):
    if(sig >= SIGMAX): sig = SIGMAX - 1.0
    bin_idx = sig / dsig
    return (int)(bin_idx)

models = ["l50n288-phew-m5", "l50n576-phew-m5"]
lgds = ["l50n288-phew-m5", "l50n576-phew-m5"]
ncpus = [256, 1024]

# models = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]
# lgds = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]
# ncpus = [128, 256]

lstyles = ["--", "-"]

# zstrs = ["z4", "z2", "z1", "z0"]
# clrs = ["yellow", "orange", "red", "purple"]
zstrs = ["z2", "z0"]
clrs = ["blue", "red"]

cols = ["atime", "ID", "mass", "siggal", "mgal", "v", "mstar", "mgas", "sfr", 'zcool', 'x', 'y', 'z']

def draw(fname, lstyle="-", clr="black"):
    df = pd.read_csv(fname, sep='\s+', names=cols)
    df['sigbin'] = pd.cut(df['siggal'], bins=sigbins, labels=sigcen)
    grps = df.groupby('sigbin')
    mtot = grps['mass'].sum()
    plt.plot(sigcen, mtot, linestyle=lstyle, color=clr)

print("Drawing for ", models[0])    
for zi, zstr in enumerate(zstrs):
    fname = DIRS['DATA']+models[0]+"/WINDS/winds."+zstr
    draw(fname, "--", clrs[zi])
print("Drawing for ", models[1])        
for zi, zstr in enumerate(zstrs):
    fname = DIRS['DATA']+models[1]+"/WINDS/winds."+zstr
    draw(fname, "-", clrs[zi])
    
plt.show()

print ("Done.")
