import pandas as pd
import numpy as np
from myinit import *
import h5py
import config
import os

def load_snapshot(snapname, grpname):
    '''
    Update the mass of existing star particles
    Create entries for new star particles
    '''
    hf = h5py.File(snapname, "r")
    header = hf['Header'].attrs
    NumPart = header['NumPart_Total']
    stars = hf['PartType4'] # Star Particles
    starIds = np.array(stars['ParticleIDs'])
    mass = np.array(stars['Masses'])

    galIds = load_galaxies(grpname, NumPart)

    sp = pd.DataFrame({'sid':starIds, 'galId':galIds, 'mass':mass})
    hf.close()

    # Need to drop duplicates
    # FUTURE: Make sure all stars have unique ID
    sp.drop_duplicates('sid', inplace=True)
    
    return sp.set_index('sid')

    # for i, sid in enumerate(starIds):
    #     if sid not in sp.index: sp.loc[sid] = {}
    #     sp.loc[sid].mass = mass[i]

def load_galaxies(fname, numPart):
    '''
    Load all SKID galaxies from a snapshot and create DataFrame gals
    Assign galId to each star particle
    '''
    i = numPart[0] + numPart[1] + 1
    galIds = []
    with open(fname, "r") as fgrp:
        while(i > 0):
            fgrp.readline()
            i -= 1
        for line in fgrp:
            galIds.append(int(line.split()[0]))
            i += 1
    assert(i == numPart[4]), "Error in load_galaxies, {} != {}.".format(i, numPart[4])
    return galIds

def update_mainId_of_stars(sp, mainIds):
    '''
    Update the MainId of all star particles in a snapshot.
      - For new stars, assign its MainId, InitMainId to its current 
    galaxy's MainId.
      - For old stars, update MainID if needed.
    '''
    mask = (sp.galId != 0)
    sp.loc[mask, 'mainId'] = mainIds[sp.loc[mask, 'galId']]
    # Leave ungrouped stars unchanged
    sp.initId.fillna(sp.mainId.astype('Int32'), inplace=True)

# sp is a pd.DataFrame that contains all stars up to this time
# sp: sid^, galId, mainId, mass
# gals: galId^, mainId, 
def find_mainId_for_gals(sp):
    '''
    Find the mainId for all SKID galaxies in a snapshot.
    Algorithm: find argmax(totMass[mainId]) for each SKID gal
      If all mainId are new, get a new mainId
    '''
    mainIds = np.zeros(sp.galId.max() + 1, dtype=int)
    for galId in sp.galId.unique():
        if galId == 0: continue
        gal = sp.query('galId == @galId')
        mtot = gal.groupby('mainId')['mass'].sum()
        if(mtot.size == 0): # new galaxy
            config.maxMainId += 1
            mainId = config.maxMainId
        else:
            mainId = mtot.index[np.argmax(mtot)]
        mainIds[galId] = mainId
    return mainIds

def process_snapshot(model, snapnum):
    snapname = os.path.join(DIRS['DATA'], model, "snapshot_{:03d}.hdf5".format(snapnum))
    grpname = os.path.join(DIRS['DATA'], model, "gal_z{:03d}.grp".format(snapnum))
    print("Load {} ...".format(snapname))
    sp = load_snapshot(snapname, grpname)
    # config.spAll['sid^', 'galId','mainId','initId'] RIGHT JOIN sp['sid^', 'mass'] ON 'sid^'
    config.spAll = pd.merge(config.spAll.drop(['galId','mass'], axis=1), sp,
                            how='right', left_index=True, right_index=True)
    print("Find mainId for galaxies ... ")
    mainIds = find_mainId_for_gals(config.spAll)
    print("Update mainId for stars... ")    
    update_mainId_of_stars(config.spAll, mainIds)
    starname = os.path.join(DIRS['DATA'], model, "stars_{:03d}.csv".format(snapnum))
    config.spAll.to_csv(starname)

def load_stars_from_snapshot(model, snapnum):
    starname = os.path.join(DIRS['DATA'], model, "stars_{:03d}.csv".format(snapnum))
    config.spAll = pd.read_csv(starname)
    config.spAll.galId = config.spAll.galId.astype('Int32')
    config.spAll.initId = config.spAll.initId.astype('Int32')
    config.spAll.mainId = config.spAll.mainId.astype('Int32')    
    return config.spAll.set_index('sid')

__mode__ = "__run__"        
if __mode__ == "__test__":
    MaxMainId = 21
    sp = pd.read_csv('example.csv')
    sp = sp.set_index('sid')
    for galId in sp.galId.unique():
        gal = sp.query('galId == @galId')
        mtot = gal.groupby('mainId')['mass'].sum()
        if(mtot.size == 0): # new galaxy
            MaxMainId += 1
            mainId = MaxMainId
        else:
            mainId = int(mtot.index[np.argmax(mtot)])
        print (galId, mainId)

if __mode__ == "__runX__":
    import matplotlib.pyplot as plt    
    config.maxMainId = 0
    config.spAll = pd.DataFrame(columns=['sid','galId','mainId','initId','mass']).set_index('sid')
    model = 'l12n144-phew-movie-200'
    # model = 'l25n288-phew-m5'
    snapnum = 25
    # config.spAll = load_stars_from_snapshot(model, snapnum)
    process_snapshot(model, snapnum)
    
    # snapname = os.path.join(DIRS['DATA'], model, "snapshot_{:03d}.hdf5".format(snapnum))
    # grpname = os.path.join(DIRS['DATA'], model, "gal_z{:03d}.grp".format(snapnum))    
    # sp = load_snapshot(snapname, grpname)
    # config.spAll = pd.merge(config.spAll.drop(['galId','mass'], axis=1), sp,
    #                  how='right', left_index=True, right_index=True)
    # mainIds = find_mainId_for_gals(config.spAll)
    # update_mainId_of_stars(config.spAll, mainIds)
    # starname = os.path.join(DIRS['DATA'], model, "stars_{:03d}.csv".format(snapnum))    

print("Compiled.")    
# note 1: All stars have been assigned a galId. For any sp.galId != 0, the galaxy must have at least one star (sp) and therefore must have a non-zero mainId.

# Largest GalId at z = 0
# 0      253.019363
# 456     69.650368
# 560     42.648006
# 101     32.455521
# 310     28.458794
