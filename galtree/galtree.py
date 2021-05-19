import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import functions
import config
from importlib import reload
import os
from myinit import *

from pygizmo import galaxy

# TODO: Generate galaxies.parquet

# For each snapshot:
# Each galId maps to a mainId
# Want to know if the mainId is
#  - consumed by any other mainId (parent)
#  - hosted by any other mainId (host)

# Comparing two stars_.csv file will tell how much mainId has changed

# Host: galId -> haloId -> hostId -> mainId(host)

def find_mainId_for_initId(sp):
    '''
    At a given snapshot, find the mapping between initId and mainId.
    initId defines where a star formed, mainId defines where it is now
    Assume most stars formed in some galaxy (initId) at earlier time end up now in a single galaxy with mainId
    '''
    initId2mainId = {}
    for initId in sp.initId.unique():
        # Get all stars that belong to the initId at this snapshot
        # i.e., all stars that originally formed in the galaxy with initId
        try:
            gal = sp.query('initId == @initId')
        except: # "initId == <NA>"
            continue
        mtot = gal.groupby('mainId')['mass'].sum()
        try:
            mainId = mtot.index[np.argmax(mtot)]
        except:
            continue
        initId2mainId[initId] = mainId
    return initId2mainId

# model = 'l25n144-phew-m5-spl'
def load_galaxy_data(path_model, snapnum, cosmo):
    '''
    Load several galaxy data of a snapshot.

    Parameters
    ----------
    path_model: str.
        Path to the data folder that contains all available data.
    cosmo: dict.
        Key, value pairs for cosmological parameters and unit conversion.

    Return
    ------
    sp: DataFrame. Star particles with their galaxy information.
    gals: DataFrame. Galaxy attributes indexed by galId
    galId2Mhost: Dict. log10(max(Mvir, Msub)) with unit conversion
    galId2hostId: Dict. Mapping from galId to hostId
    '''
    parname = os.path.join(path_model, "so_z{:03d}.par".format(snapnum))
    host = pd.read_csv(parname, sep='\s+', names=['galId', 'hostId']).set_index('galId')
    galId2hostId = host.to_dict()['hostId']
    path_stat = os.path.join(path_model, "gal_z{:03d}.stat".format(snapnum))
    path_sovcirc = os.path.join(path_model, "so_z{:03d}.sovcirc".format(snapnum))
    gals = galaxy.read_stat(path_stat)
    halos = galaxy.read_sovcirc(path_sovcirc)
    galId2Mhost = halos.apply(lambda x: np.log10(np.maximum(x.Mvir, x.Msub)/cosmo['h']), axis=1).to_dict()
    galId2Mgal = gals[['Mtot', 'Mstar']].apply(lambda x: np.log10(x * cosmo['msun_tipsy'] / cosmo['h']))
    return galId2Mgal, galId2Mhost, galId2hostId

def find_galId_for_mainId(sp):
    '''
    Find the galId that corresponds to each mainId at a snapshot.
    One mainId may relate to multiple galId, so we pick the most massive one.
    
    Parameters
    ----------
    sp: DataFrame
        All stars at this snapshot, each having a galId and a mainId.

    Return
    ------
    mainId2galId
    '''
    # Only select stars that are in a galaxy (galId > 0)
    sp = sp[sp.galId > 0]
    # Group by mainId, and pick the galId with the largest sum(mass)
    grps = sp[['mainId','galId','mass']].groupby(['mainId','galId']).sum('mass')
    grps = grps.reset_index()
    idx = grps.groupby('mainId')['mass'].transform(max) == grps['mass']
    mainId2galId = grps[idx][['mainId', 'galId']].set_index('mainId')
    return mainId2galId

def build_mainId_map_for_snapshot(path_model, snapnum):
    galId2Mgal, galId2Mhost, galId2hostId = load_galaxy_data(path_model, snapnum, cosmo)
    sp = functions.load_stars_from_snapshot(path_model, snapnum)

    mainIdMap = find_galId_for_mainId(sp)
    # Get Mtot, Mstar, Mhost for the table
    mainIdMap['hostId'] = mainIdMap['galId'].map(galId2hostId)
    mainIdMap = pd.merge(mainIdMap, galId2Mgal, how='left', left_on='galId', right_index=True)
    mainIdMap['Mhost'] = mainIdMap['hostId'].map(galId2Mhost)
    mainIdMap['snapnum'] = snapnum
    return mainIdMap

def find_main_descendent_galId(sp, spn):
    '''
    For each galaxy from a snapshot, find its main descendent in the next snapshot. The main descendent of a galaxy is defined as the galaxy that has a majority of its stellar mass.

    Parameters
    ----------
    sp: DataFrame.
        Star particles from this snapshot.
    spn: DataFrame.
        Star particles from the next snapshot.

    Returns
    -------
    galId2galIdNext: dict.
        A temporary mapping between galId to galIdNext.
    '''
    sppair = pd.merge(sp[sp.galId > 0][['galId','mass']],
                      spn[['galId']].rename({'galId':'galIdNext'}, axis=1),
                      how='left', left_index=True, right_index=True)
    grp = sppair.groupby(['galId', 'galIdNext']).sum('mass').reset_index()
    idx = grp.groupby('galId')['mass'].transform(max) == grp['mass']
    return dict(zip(grp[idx].galId, grp[idx].galIdNext))

model = 'l12n144-phew-movie-200'
path_model = os.path.join(DIRS['DATA'], model)
cosmo = {'h':0.70, 'msun_tipsy':4.3e16}
mainIdMap = build_mainId_map_for_snapshot(path_model, 100)

sp = functions.load_stars_from_snapshot(path_model, 100)
spn = functions.load_stars_from_snapshot(path_model, 105)
galId2galIdNext = find_main_descendent_galId(sp, spn)
mainIdMap['mainIdNext'] = mainIdMap['galId'].map(galId2galIdNext)

# How to define parent

# construct mainId table:
# ISSUE: many galId -> same mainId
# SOLUTION: Choose the largest gal
# sp[['galId', 'mainId', 'mass']]
# groupby 'mainId', 'galId'
# agg: sum(mass)
# pick galId with largest sum(mass)

# show_merger_tree(snapnum, galId)
# 1. query(snapnum, galId) from the stars_{}.csv file
#    -> List of stars in the galaxy
# 2. GroupBy initId of all particles in it.
# 3. For each initId, plot the mtot of mainId = initid vs. snapnum
#    Ending with the parent's mtot



# [Code Recycle Bin]
# initId2mainId = find_mainId_for_initId(sp)

# df=pd.DataFrame({'InitId':initId2mainId.keys(),'mainId':initId2mainId.values()})

# gals = sp[['galId','mainId']].drop_duplicates('galId').set_index('galId')
# galId2mainId = gals.to_dict()['mainId']
# host['mainIdHost'] = host.hostId.map(galId2mainId).astype('Int32')
# gals = pd.merge(gals, host.drop('hostId',axis=1), how='left', left_index=True, right_index=True)
# gals.mainIdHost.fillna(gals.mainId, inplace=True)

