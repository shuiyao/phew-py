import h5py
from myinit import *
import numpy as np
import pandas as pd
import time

# Find the host haloes for PhEW particles.
# Input: initwinds.*, .sogrp, particleIDs and masses
# Output: idx PhEWKey HID Msub

tbeg = time.time()

cols_halos = ["HID", "Mvir", "Rvir", "R_vmax", "Vmax", "Npart", "Msub", "Rsub"]
dtypes_halos = {'HID':'int', 'Npart':'int'}

model = "l50n576-phew-m5"
AMAX = 0.505

snapstr = "078"
path_model = os.path.join(DIRS['DATA'], model)
path_hdf5 = os.path.join(path_model, "snapshot_{}.hdf5".format(snapstr))
path_sogrp = os.path.join(path_model, "so_z{}.sogrp".format(snapstr))
path_sovcirc = os.path.join(path_model, "so_z{}.sovcirc".format(snapstr))
path_output = os.path.join(path_model, "snapshot_{}.phewsHalos".format(snapstr))

path_schema = os.path.join(path_model, "WINDS/initwinds.0")
path_initwinds = os.path.join(path_model, "WINDS/initwinds.")

def read_hdf5():
    hf = h5py.File(path_hdf5, "r")
    gp = hf['PartType0'] # Gas Particles
    n_gas = hf['Header'].attrs['NumPart_Total'][0]
    pdDF = pd.DataFrame({"idx":range(n_gas), "PID":gp['ParticleIDs'], 'MassHDF5':gp['Masses']})
    return pdDF

def read_haloes():
    df = pd.read_csv(path_sovcirc, sep='\s+', skiprows=1,
                     names=cols_halos, dtype=dtypes_halos)
    gids = pd.read_csv(path_sogrp, sep='\s+', header=0, names=['HID'])
    df = df[:-1].set_index("HID")
    return df, gids

def read_initwinds(ncpu=256):
    df = pd.read_csv(path_schema, sep='\s+', header=0)
    schema = df.columns
    df = df[['#atime', 'ID', 'Mass', 'PID']]
    for i in range(ncpu)[1:]:
        dfnew = pd.read_csv(path_initwinds+str(i), sep='\s+', names=schema)
        dfnew = dfnew[['#atime', 'ID', 'Mass', 'PID']]
        df = pd.concat([df, dfnew])
    return df

#raise ValueError()

print("Loading data...")
df_halos, gids = read_haloes()

df = read_hdf5()
df_winds = read_initwinds()
df_winds.rename(columns={'#atime':"atime", "ID":"PhEWKey"}, inplace=True)
df_winds = df_winds.query('atime < @AMAX and atime > 0.50')

print("Matching PhEW to gas particles...")
df2 = df_winds.merge(df, left_on='PID', right_on='PID', how='inner')

# Release Memory
del df
del df_winds

df2['MassDiff'] = abs(df2['Mass'] - df2['MassHDF5'])
df2 = df2.groupby('PID').min('MassDiff')

print("Matching Haloes...")
df2 = df2[['idx', 'PhEWKey']].merge(gids, how='left', left_on='idx', right_index=True)
df2 = df2.query("HID > 0")
df2 = df2.merge(df_halos[['Mvir', 'Rvir', 'Msub']], how='left', left_on='HID', right_index=True)
df2['LogMvir'] = np.log10(df2['Mvir'] / 0.7)
df2['LogMsub'] = np.log10(df2['Msub'] / 0.7)
df2.drop(['Mvir', 'Msub'], axis=1).to_csv(path_output)

# fout.write("#idx ID HID LogMvir LogMsub Rvir\n")

print("Cost: {}".format(time.time() - tbeg))
print('DONE')

# Cost: 19.532994747161865
