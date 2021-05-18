# Input: snapshot_*.sph
# - bash loadhdf5.sh $modelname $redshift:
#   - get_particles -sph $modelname
# - Use NOSPHSKIP to include every single SPH particle in it

import pandas as pd
import numpy as np
from myinit import *
import seaborn as sns

import matplotlib.pyplot as plt    

mode = "match"
Subset = "all"

#model = "l50n288-phew-m5"
model = "l50n288-phewoff"
zstr = "108"
folder = os.path.join(DIRS['DATA'], model)
Fgas = os.path.join(folder, "snapshot_"+zstr+".sph")
Fsogrp = os.path.join(folder, "so_z"+zstr+".sogrp")
Fsovirc = os.path.join(folder, "so_z"+zstr+".sovcirc")
Fout = model+"_"+zstr+".Zhist"

ZSOLAR = np.log10(0.0122)    
NROWS = None

def write_Zhist_to_csv(hist, subset='all'):
    x, y = hist.get_lines()[0].get_data()
    dfout = pd.DataFrame({"OH":x})
    i = 1
    # BE CAREFUL ABOUT THE WAY lines ARE ORDERED.
    for line in hist.get_lines()[:]:
        _, y = line.get_data()
        dfout['mbin'+str(i)] = y
        i = i + 1
    if(subset=='cold'):
        dfout.to_csv(Fout+"_cold", index=False)
    elif(subset=='hot'):
        dfout.to_csv(Fout+"_hot", index=False)
    elif(subset=='all'):
        dfout.to_csv(Fout, index=False)        
        

if(mode == "load"):
    df = pd.read_csv(Fgas, sep='\s+', nrows=NROWS)
    h = pd.read_csv(Fsovirc, sep='\s+')
    h = h[['#', 'grp#', 'Npart']]
    h = h.rename(columns={'#':"galId", 'grp#':'Mvir', 'Npart':'Msub'}).set_index('galId')
    gids = pd.read_csv(Fsogrp, sep='\s+', nrows=NROWS, header=0)
    df['galId'] = gids
    df = df[df.galId > 0]
    df = df[df.SfFlag == 0]
    print("Loaded.")

if(mode == "match"):
    cen = h[h.Msub > 0]
    cen['logM'] = np.log10(cen['Msub'])
    hZ = df[['galId','logT','Mass','MassZ']].join(cen[['logM']], how='inner', on='galId')
    hZ['OH'] = np.log10(hZ['MassZ']) - ZSOLAR
    # bins = np.linspace(10., 13.5, 41)
    bins = np.linspace(10.5, 14.5, 5)
    bcen = 0.5*(bins[1:] + bins[:-1])
    mbins = pd.cut(hZ.logM, bins=bins, labels=bcen)
    hZ['Mbin'] = mbins
    HotIdx = hZ['logT'] > 5.5
    Hot = hZ[HotIdx]
    Cold = hZ[~HotIdx]
    print("Plotting Histogram...")
    if(Subset == 'cold'):
        hZ2 = Cold[Cold['OH'] > -6.0]
    elif(Subset == 'hot'):        
        hZ2 = Hot[Hot['OH'] > -6.0]
    elif(Subset == 'all'):
        hZ2 = hZ[hZ['OH'] > -6.0]
    # hist = sns.kdeplot(data=hZ2, x='OH', hue='Mbin', palette='rocket', common_norm=False)
    # hist.get_lines[:].get_data()
    Zbins = linspace(-4.0, 1.0, 26)
    hist = sns.histplot(data=hZ2, x='OH', bins=Zbins, hue='Mbin', element='step',fill=False, palette='rocket', hue_order=[11.0,12.0,13.0,14.0])
    write_Zhist_to_csv(hist, subset=Subset)
    plt.close()

