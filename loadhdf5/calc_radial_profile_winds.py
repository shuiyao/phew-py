# Calculate the cumulative distribution of wind mass as a function of radius
# Divide the wind into old and younger winds

from myinit import *
import pandas as pd
from cosmology import tcosmic, acosmic
import scipy
import numpy as np
import seaborn as sns

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 10
nagebins = 5
nsigbins = 5
redges = linspace(0., 1., nbins+1)
dr = 1. / nbins

#method = 'fill'
multiple_method = 'stack'

bins_rad = linspace(0.0, 1.0, nbins+1)
bins_rad_label = 0.5 * bins_rad[1:] + 0.5 * bins_rad[:-1]
bins_age = linspace(0.0, 10.0, nagebins+1)
bins_age_label = 0.5 * bins_age[1:] + 0.5 * bins_age[:-1]
bins_sig = linspace(0.0, 150.0, nsigbins+1)
bins_sig_label = 0.5 * bins_sig[1:] + 0.5 * bins_sig[:-1]

def find_radial_bin(r):
    if(r > 1.0): return -1
    bin_idx = r / dr
    return (int)(bin_idx)

REDSHIFT = 0.25
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 0.25): zstr = "098"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"

tcurrent = tcosmic(1./(REDSHIFT+1.0))

TMAXCUT = 5.0

model = 'l50n576-phew-m5'
# model = 'l50n288-phew-m5'

def calculate_gas_profile(model, zstr, mstr, PhEW=True):

    fname = DIRS['SCIDATA'] + model + "/snapshot_"+zstr+".gas."+mstr
    fauxname = DIRS['SCIDATA'] + model + "/snapshot_"+zstr+".gasaux."+mstr        
    gas = pd.read_csv(fname, sep='\s+', header=0)
    gasaux = pd.read_csv(fauxname, sep='\s+', header=0)
    # gas.rename(columns={'#Idx':'Idx'})
    # gasaux.rename(columns={'#f_w':'f_w'})
    gas = gas[gas['SfFlag'] == 0]
    concat = pd.concat([gas[['WMass', 'dr', 'Rvir']], gasaux[['t', 'sig']]], axis=1)
    # concat['rbin'] = concat.apply(lambda x : find_radial_bin(x.dr / x.Rvir), axis=1)
    concat['rbin'] = concat['dr'] / concat['Rvir']

    concat['t'] = (tcurrent - concat['t']) / 1.e9
    concat.rename(columns={'t':'age'}, inplace=True)

    rbins = pd.cut(concat['rbin'], bins=bins_rad, labels=bins_rad_label)
    agebins = pd.cut(concat['age'], bins=bins_age, labels=bins_age_label)
    sigbins = pd.cut(concat['sig'], bins=bins_sig, labels=bins_sig_label)    
    bin_sum = concat.groupby([rbins, agebins]).sum()

    # concat['weighted_age'] = concat['WMass'] * concat['age']
    # bin_sum = concat[['WMass', 'weighted_age']].groupby(rbins).sum()
    # bin_sum['weighted_age'] = tcurrent - bin_sum['weighted_age'] / bin_sum['WMass']

    print ("Reading: ", fname)
    return bin_sum, concat

import matplotlib.pyplot as plt

# binsum, concat = calculate_gas_profile("l50n288-phew-m5", zstr, 'mh13')

# calculate_gas_profile("l50n288-phew-m4", zstr)
# calculate_gas_profile("l25n288-phew-m5", zstr)
#calculate_gas_profile("l50n288-phewoff", zstr, PhEW=False)
#binsum = calculate_gas_profile("l50n576-phew-m5", zstr)

# ncells_x = nbins # rbins
# ncells_y = nagebins
# ncells = ncells_x * ncells_y
# xnodes = bins_rad_label
# ynodes = bins_age_label
# z = array(binsum['WMass'])
# z = np.log(z)
# x, y = scipy.meshgrid(xnodes, ynodes)
# z = scipy.reshape(z, (ncells_x, ncells_y)).T

# fig, axs = plt.subplots(2,1,figsize=(8,8))
# axs.flatten()
# axs[0].pcolor(xnodes, ynodes, z, cmap=plt.cm.Purples)
sns.set_theme(style='ticks')
fig, axs = plt.subplots(3,1, figsize=(7,9))
axs = axs.flatten()
titles = [r"$11.0 < \log(M_{vir}/M_\odot) < 11.5$",\
          r"$11.85 < \log(M_{vir}/M_\odot) < 12.15$",\
          r"$12.85 < \log(M_{vir}/M_\odot) < 13.15$",]

feature, hue = "age", "Age [Gyr]"
bins, bins_label = bins_age, bins_age_label

# feature, hue = "sig", "sigma [km/s]"
# bins, bins_label = bins_sig, bins_sig_label

for i, mstr in enumerate(["mh11", "mh12", "mh13"]):
    ax = axs[i]
    binsum, concat = calculate_gas_profile(model, zstr, mstr)
    # ybins = pd.cut(concat[feature], bins=bins_age, labels=bins_age_label)
    ybins = pd.cut(concat[feature], bins=bins, labels=bins_label)    
    concat[feature] = ybins
    concat.rename(columns={feature:hue}, inplace=True)
    lgd = True if i == 0 else False

    print("Making Plot...")
    sns.despine(fig)
    sns.histplot(concat, x="rbin", bins=bins_rad, weights='WMass', hue=hue, \
                 multiple=multiple_method, palette="rocket", \
                 edgecolor=".3", linewidth=.5, stat='density',\
                 hue_order = bins_label[::-1], ax=ax, legend=lgd)
    # ax.set_ylabel("Fraction")
    ax.set_ylabel("Mwind")
    ax.set_title(titles[i])
    ax.set_xlim(0.0, 1.0)
    if(i != 2):
        ax.set_xticks([])
        ax.set_xlabel("")
    
axs[2].set_xlabel(r"$r/R_{vir}$")
plt.savefig(DIRS['FIGURE']+"tmp.pdf")

print ("Done.")
