# Input: snapshot_*.sph
# - bash loadhdf5.sh $modelname $redshift:
#   - get_particles -sph $modelname
# - Use NOSPHSKIP to include every single SPH particle in it

import pandas as pd
import numpy as np
from myinit import *
import seaborn as sns

import matplotlib.pyplot as plt    

mode = "save"
# testplot, plot

#model = "l50n288-phew-m5"
models = ["l50n288-phew-m5", "l50n288-phewoff", "l50n576-phew-m5"]

zstr = "098"

ZSOLAR = np.log10(0.0122)    
NROWS = None

def load_model(model, zstr, mstr):
    Fgas = os.path.join(scifolder, "snapshot_"+zstr+".gas."+mstr)
    Fsovirc = os.path.join(folder, "so_z"+zstr+".sovcirc")
    print('Reading: ', Fgas)
    df = pd.read_csv(Fgas, sep='\s+')
    h = pd.read_csv(Fsovirc, sep='\s+')
    h = h[['#', 'grp#', 'Npart']]
    h = h.rename(columns={'#':"galId", 'grp#':'Mvir', 'Npart':'Msub'}).set_index('galId')
    # Central halo only
    h = h[h.Msub > 0]
    h['logM'] = np.log10(h['Msub'])
    # Clean the gas particle data set
    print("Total number of particles: {:6d}".format(df.shape[0]))
    df = df.query('HID > 0')
    print("In Halos: {:6d}".format(df.shape[0]))
    df = df.query('SfFlag == 0')
    print("Non Star-Forming: {:6d}".format(df.shape[0]))    
    df = df.query('dr > 0.1 * Rvir')
    print("r > 0.1 Rvir: {:6d}".format(df.shape[0]))
    hZ = df[['HID','logT','Mass','Z']].join(h[['logM']], how='inner', on='HID')
    print("In central haloes: {:6d}".format(hZ.shape[0]))
    # hZ = hZ.query('Z > 0')
    hZ['OH'] = np.log10(hZ['Z']+1.e-9) - ZSOLAR
    # print("Non-zero metallicty: {:6d}".format(hZ.shape[0]))
    return hZ


def show_model(hZ, ax=None):
    hZ['Hot'] = (hZ['logT'] > 5.0)
    Zbins = linspace(-4.0, 1.0, 51)
    hist = sns.histplot(data=hZ, x='OH', bins=Zbins, hue='Hot', element='step', fill=False, ax=ax)

if(mode == "test"):
    # Z < 1.e-6: phew-m5: 0.57%
    # Z < 1.e-6: phewoff: 13.6%
    
    fig, axs = plt.subplots(3,1,figsize=(6,9))
    mstrs = ['mh11', 'mh12', 'mh13']
    for i, mstr in enumerate(mstrs):
        folder = os.path.join(DIRS['DATA'], model)        
        scifolder = os.path.join(DIRS['SCIDATA'], model)    
        hZ = load_model(model, zstr=zstr, mstr=mstr)
        show_model(hZ, ax=axs[i])
    plt.savefig(DIRS["FIGURE"]+"tmp.pdf")

def get_metal_histogram(hZ, zbins=50):
    Zbins = linspace(-4.0, 1.0, 51)    
    hZ['Hot'] = (hZ['logT'] > 5.0)
    Zhot = hZ.query('Hot==True')[['OH','Mass']]
    Zcold = hZ.query('Hot==False')[['OH','Mass']]
    yhot, _ = np.histogram(Zhot.OH, bins=Zbins, weights=Zhot.Mass)
    ycold, _  = np.histogram(Zcold.OH, bins=Zbins, weights=Zcold.Mass)
    yhot = np.insert(yhot, 0, Zhot.query('OH <= -4.0').Mass.sum())
    ycold = np.insert(ycold, 0, Zcold.query('OH <= -4.0').Mass.sum())
    ytot = yhot.sum() + ycold.sum()
    yhot /= ytot
    ycold /= ytot
    return ycold, yhot

def save_model(model):
    # For each model, want cold/hot, mh11, mh12, mh13
    Fout = os.path.join(scifolder, "Zhist_"+zstr)
    Zbins = linspace(-4.0, 1.0, 51)
    Zcen = 0.5 * (Zbins[1:] + Zbins[:-1])
    Zcen = np.insert(Zcen, 0, -4.0)
    dfout = pd.DataFrame({'logZ':Zcen})
    mstrs = ['mh11', 'mh12', 'mh13']
    for i, mstr in enumerate(mstrs):
        hZ = load_model(model, zstr=zstr, mstr=mstr)
        c, h = get_metal_histogram(hZ, zbins=Zbins)
        dfout[mstr+'c'] = c
        dfout[mstr+'h'] = h
    print("Writing File: ", Fout)
    dfout.to_csv(Fout, index=False)

if(mode == "save"):
    for model in models:
        folder = os.path.join(DIRS['DATA'], model)
        scifolder = os.path.join(DIRS['SCIDATA'], model)
        save_model(model)
