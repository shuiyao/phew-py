from mymod import *
import matplotlib.pyplot as plt
import config_mpl
from cosmology import acosmic, tcosmic
from numpy import isnan, histogram, cumsum

Mgas_lres = 0.00653561
Mgas_hres = 0.00081695

# models = ["l25n144-phew-m5-spl"]
# lgds = ["l25n144-spl"]
models = ["l25n144-phew-m5"]
lgds = ["l25n144"]
ls = ["-"]

modeli = 0
mstr = "mh11"
zstr = "108"
folder = "/scratch/shuiyao/scidata/gadget3io/"
fname = folder + models[modeli] + "_"+zstr+".gas."+mstr
fnameaux = folder + models[modeli]+"/"+models[modeli]+"_"+zstr+".gasaux." + mstr

def get_color(vals, vmin=0.0, vmax=1.0, cmap=plt.get_cmap('jet')):
    cgrad = vmax - vmin
    clrsarr = []
    for val in vals:
        cval = (val - vmin) / cgrad
        if(val < vmin): cval = 0.0
        if(val > vmax): cval = 1.0
        clrsarr.append(cmap(cval))
    return clrsarr

def reverse_cumsum(hist, edges):
    cen = 0.5 * (edges[1:] + edges[:-1])
    cumfrac = cumsum(hist[::-1])[::-1]
    cumfrac = cumfrac / (float)(cumfrac[0])
    return cen, cumfrac

def show_histogram():
    fig, ax = plt.subplots(1,1, figsize=(8,7))
    gasaux = genfromtxt(fnameaux, names=True)
    # gas = genfromtxt(fname, names=True)
    # awind = []
    # for i in range(len(gasaux)):
    #     if(gasaux['t'][i] < 1.34e10):
    #         awind.append(acosmic(gasaux['t'][i]))
    #     else: awind.append(1.0)
    # hist, edges = histogram(awind, bins=100)
    gasaux = gasaux[gasaux['sig'] < 150.]
    hist, edges = histogram(gasaux['sig'], bins=linspace(0., 150., 100))
    cen, cumfrac = reverse_cumsum(hist, edges)
    ax.plot(cen, cumfrac, "b", linestyle=ls[modeli])
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")    
    plt.show()

def show_scatter_plot():
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    gas = genfromtxt(fnameaux, names=True)
    gas = gas[gas['var_t'] > 0]
    gas = gas[gas['t']+gas['var_t']<1.34e10]
    gas = gas[gas['t']-gas['var_t']>1.e6]
    gas = gas[::100]
    clrs = get_color(gas['f_w'])

    # atime = acosmic(gas['t'])
    # err1_a = atime - acosmic(gas['t'] - gas['var_t'])
    # err2_a = acosmic(gas['t'] + gas['var_t']) - atime

    cmap = plt.get_cmap('jet')
    # ax.errorbar(atime, gas['sig'], xerr=[err1_a, err2_a], yerr=gas['var_sig'], ls='none', color='black')
    # ax.scatter(atime, gas['sig'], s=6, c=clrs, cmap=cmap)

    ax.errorbar(gas['t']/1.e9, gas['sig'], xerr=gas['var_t']/1.e9, ls='none', color='black')    
    # ax.errorbar(gas['t']/1.e9, gas['sig'], yerr=gas['var_sig'], ls='none', color='black')
    ax.scatter(gas['t']/1.e9, gas['sig'], s=12, c=clrs, cmap=cmap)
    ax.set_xlabel(r"$t_{wind} [Gyr]$")
    ax.set_ylabel(r"$\sigma_{gal} [km/s]$")
    ax.set_title(models[modeli]+", Low-Mass")
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()

print "done."
