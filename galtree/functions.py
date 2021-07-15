import pandas as pd
import numpy as np
from myinit import *
import h5py
import config
import os

schema_columns = ['sid','galId','mainId', 'initId', 'mass']
schema_dtypes = {'sid': "int64",
                 'galId': "int32",
                 'mainId': "int32",
                 'initId': "int32",
                 'mass': "float32"}

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

def reset_mainId_for_rogue_and_traitor_stars(sp, galIdPrev2galId):
    '''
    Set the mainID for all rogue stars (galId = 0) and traitor stars (not belonging to the main descendent of the mainId).

    galIdPrev2galId: Obtained from find_main_descendent_galId(sp)
    '''
    # Traitor star:  sp.galId != galId2galIdNext[sp.galIdPrev]
    # Rogue star: sp.galId == 0.
    mask = (sp.galId != 0) & (sp.galId == sp.galIdPrev.map(galIdPrev2galId))
    sp.loc[~mask, 'mainId'] = 0

def update_mainId_and_initId_of_stars(sp, mainIds):
    '''
    Update the MainId of all star particles in a snapshot.
      - For new stars, assign its MainId, InitMainId to its current 
    galaxy's MainId.
      - For old stars, update MainID if needed.
    '''
    mask = (sp.galId != 0) # Leave ungrouped stars unchanged
    sp.loc[mask, 'mainId'] = mainIds[sp.loc[mask, 'galId']]
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
    # mainIds = {}
    for galId in sp.galId.unique(): # cover all galaxies having stars
        if galId == 0: continue
        galsp = sp[sp.mainId > 0].query('galId == @galId')
        mtot = galsp.groupby('mainId')['mass'].sum()
        if(mtot.size == 0): # new galaxy
            config.maxMainId += 1
            mainId = config.maxMainId
        else:
            mainId = mtot.index[np.argmax(mtot)]
        mainIds[galId] = mainId
    return mainIds

def process_snapshot(model, snapnum, output=True):
    snapname = os.path.join(DIRS['DATA'], model, "snapshot_{:03d}.hdf5".format(snapnum))
    grpname = os.path.join(DIRS['DATA'], model, "gal_z{:03d}.grp".format(snapnum))
    print("Load {} ...".format(snapname))
    sp = load_snapshot(snapname, grpname)

    # Update galIdPrev for future use (identify main descendent)
    config.spAll['galIdPrev'] = config.spAll['galId']
    # Drop galId and Mass for existing stars and add new stars to the global spAll
    # Essentially: config.spAll['sid^', 'galIdPrev','mainId','initId']
    #                          RIGHT JOIN sp['sid^', 'mass'] ON 'sid^'
    config.spAll = pd.merge(config.spAll.drop(['galId','mass'], axis=1), sp,
                            how='right', left_index=True, right_index=True)

    galIdPrev2galId = find_main_descendent_galId(config.spAll)

    reset_mainId_for_rogue_and_traitor_stars(config.spAll, galIdPrev2galId)
    
    print("Find mainId for galaxies ... ")
    mainIds = find_mainId_for_gals(config.spAll)

    print("Update mainId for stars... ")    
    update_mainId_and_initId_of_stars(config.spAll, mainIds)

    if(output == True):
        starname = os.path.join(DIRS['DATA'], model, "stars_{:03d}.csv".format(snapnum))
        config.spAll.drop('galIdPrev', axis=1).to_csv(starname)

def find_main_descendent_galId(sp):
    '''
    For each galaxy from a snapshot, find its main descendent in the next snapshot. The main descendent of a galaxy is defined as the galaxy that has a majority of its stellar mass.
    Algorithm:
    1. Select stars that was found in a galaxy in both the previous and the current snapshot (galIdPrev > 0 and galId > 0)
       - There might be galIdPrev that do not have any corresponding galId. These galaxies are lost (tidally destroyed) in the current snapshot.
    2. Find the total mass of stars, mtot, in any galId that was from galIdPrev
    3. Pick the galId with max(mtot) as the descendent of galIdPrev

    Parameters
    ----------
    sp: DataFrame.
        Star particles from this snapshot.
        Should have ['galIdPrev', 'galId', 'mass']

    Returns
    -------
    galIdPrev2galId: dict.
        A temporary mapping between galIdPrev to galId.
        Note: not all galId2Prev from sp has a mapping.
    '''
    # sppair = pd.merge(sp[sp.galId > 0][['galId','mass']],
    #                   spn[['galId']].rename({'galId':'galIdNext'}, axis=1),
    #                   how='left', left_index=True, right_index=True)

    sp = sp[(sp.galIdPrev > 0) & (sp.galId > 0)]
    grp = sp.groupby(['galIdPrev', 'galId']).sum('mass').reset_index()
    # >>> (galIdPrev, galId), sum(mass)
    idx = grp.groupby('galIdPrev')['mass'].transform(max) == grp['mass']
    # 
    return dict(zip(grp[idx].galIdPrev, grp[idx].galId))
    # tmp = spn.drop_duplicates('galId')
    # galId2mainId = dict(zip(tmp.galId, tmp.mainId))
    # galId2mainId[0] = 0 # Rogue stars
    # return dict(zip(grp[idx].galId, grp[idx].galIdNext.map(galId2mainId)))

def load_stars_from_snapshot(path_model, snapnum):
    starname = os.path.join(path_model, "stars_{:03d}.csv".format(snapnum))

    # config.spAll = pd.read_csv(starname, names=schema_columns, dtype=schema_dtypes)
    config.spAll = pd.read_csv(starname)
    config.spAll.galId = config.spAll.galId.astype('Int32')
    config.spAll.initId = config.spAll.initId.astype('Int32')
    config.spAll.mainId = config.spAll.mainId.astype('Int32')    
    return config.spAll.set_index('sid')

__mode__ = "__runX__"        
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

if __mode__ == "__run__":
    import matplotlib.pyplot as plt    
    config.maxMainId = 0
    config.spAll = pd.DataFrame(columns=schema_columns).set_index('sid')
    model = 'l12n144-phew-movie-200'
    # model = 'l25n288-phew-m5'

    path_model = os.path.join(DIRS['DATA'], model)
    config.spAll = load_stars_from_snapshot(path_model, 198)
    print("processing new snapshot.")
    process_snapshot(model, 200, output=False)
    
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
