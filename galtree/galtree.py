import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import functions
import config
from importlib import reload
import os
from myinit import *
import glob

from pygizmo import galaxy

# TODO: Generate galaxies.parquet

# For each snapshot:
# Each galId maps to a mainId
# Want to know if the mainId is
#  - consumed by any other mainId (parent)
#  - hosted by any other mainId (host)

# Comparing two stars_.csv file will tell how much mainId has changed

# Host: galId -> haloId -> hostId -> mainId(host)

schema_columns = ['mainId', 'snapnum', 'galId', 'hostId', 'Mtot', 'Mstar', 'Mhost', 'mainIdNext']
schema_dtypes = {'mainId': "int32",
                 'snapnum': "int32",
                 'galId': "int32",
                 'hostId': "int32",
                 'Mtot': "float32",
                 'Mstar': "float32",
                 'Mhost': "float32",                 
                 'mainIdNext': "int32"}

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
    galId2Mgal = gals[['Mtot', 'Mstar']].apply(lambda x: np.log10(x * cosmo['msun_tipsy'] / cosmo['h']))
    halos = galaxy.read_sovcirc(path_sovcirc)
    galId2Mhost = halos.apply(lambda x: np.log10(np.maximum(x.Mvir, x.Msub)/cosmo['h']), axis=1).to_dict()
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

def find_main_descendent_mainId(sp, spn):
    '''
    For each galaxy (galId, mainId) in the current snapshot, find the mainId^ of its main descendent galaxy in the next snapshot. First step is to find the galId^ of the main descendent. There are several cases:
    I. There is no galId^ to be found: 
       The galaxy has been tidally distroyed. Set mainIdNext to 0. 
    II. Main descendent exists, but most of its mass came from another galaxy.
       The galaxy has merged into a parent galaxy. Set mainIdNext to mainId^
    III. Main descendent exists, and it is the same galaxy.
       Simply set mainIdNext to mainId^, which is equal to mainId.
    
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
        Note: not all galId2galIdNext from sp has a mapping.
    '''
    sp = pd.merge(sp[sp.galId > 0][['galId','mass']],
                  spn[spn.galId > 0][['galId']].rename({'galId':'galIdNext'}, axis=1),
                  how='left', left_index=True, right_index=True)

    grp = sp.groupby(['galId', 'galIdNext']).sum('mass').reset_index()
    idx = grp.groupby('galId')['mass'].transform(max) == grp['mass']

    tmp = spn.drop_duplicates('galId')
    galId2mainId = dict(zip(tmp.galId, tmp.mainId))
    galId2mainId[0] = 0 # Rogue stars
    return dict(zip(grp[idx].galId, grp[idx].galIdNext.map(galId2mainId)))

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

model = 'l12n144-phew-movie-200'
path_model = os.path.join(DIRS['DATA'], model)
cosmo = {'h':0.70, 'msun_tipsy':4.3e16 * (12./25.)**3}
list_of_snapnums = list(range(0, 201, 2))

glob.glob(path_model+"/stars_*.csv")

def build_mainId_table(path_model, cosmo, list_of_snapnums):
    for i, snapnum in enumerate(list_of_snapnums):
        print("snapnum = {}".format(snapnum))
        galId2Mgal, galId2Mhost, galId2hostId = load_galaxy_data(path_model, snapnum, cosmo)
        if(i == 0):
            sp = functions.load_stars_from_snapshot(path_model, snapnum)
        else:
            sp = spn.copy()
        if(snapnum != list_of_snapnums[-1]): # not the last one
            spn = functions.load_stars_from_snapshot(path_model, list_of_snapnums[i+1])
            
        mainIdMap = find_galId_for_mainId(sp)
        # Get Mtot, Mstar, Mhost for the table
        mainIdMap['hostId'] = mainIdMap['galId'].map(galId2hostId)
        mainIdMap = pd.merge(mainIdMap, galId2Mgal,
                             how='left', left_on='galId', right_index=True)
        mainIdMap['Mhost'] = mainIdMap['hostId'].map(galId2Mhost)
        mainIdMap['snapnum'] = snapnum
        if(snapnum != list_of_snapnums[-1]):
            galId2mainIdNext = find_main_descendent_mainId(sp, spn)
            # Not every galId can be mapped.
            # The galIds not in the dictionary map to NaN. As a result, the field mainIdNext is automatically treated as float64
            mainIdMap['mainIdNext'] = mainIdMap['galId'].map(galId2mainIdNext).astype('Int32').fillna(0)
        else:
            mainIdMap['mainIdNext'] = 0
        mainIdTable = mainIdMap.copy() if (i == 0) else pd.concat([mainIdTable, mainIdMap])
    return mainIdTable

# mainIdTable = build_mainId_table(path_model, cosmo, list_of_snapnums)
# mainIdTable.reset_index().to_csv('test.csv', columns=schema_columns, index=False)

x = pd.read_csv('test.csv', header=0, dtype=schema_dtypes)
sp = functions.load_stars_from_snapshot(path_model, 200)
gals = x[x.snapnum==200]

def code_repo():
    # Get the mainIds of the most massive galaxies at z = 0
    bgals = x[x.snapnum==200].sort_values('Mstar', ascending=False).iloc[:5]
    mgals = gals[(10.8 < gals.Mstar) & (gals.Mstar < 11.2)]

    import matplotlib.pyplot as plt
    # Let's try with the most massive galaxy: 4813
    mainIdTarget = 20679

    gal = x[x.mainId == mainIdTarget]
    plt.plot(gal.snapnum, gal.Mstar, "b-")

    upstr = x[x.mainIdNext == mainIdTarget].drop_duplicates('mainId') # 934 rows
    upstr = upstr[(upstr.mainId != mainIdTarget) & (upstr.Mstar > 10.0)] # 50 rows

# sub = gals[gals.Mstar > 10.0]
# for mainId in sub.mainId:
#     gal = x[x.mainId == mainId]
#     plt.plot(gal.snapnum, gal.Mtot, "r-", alpha=0.5)
#     plt.plot(gal.snapnum, gal.Mhost, "b-", alpha=0.5)
# plt.show()














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

