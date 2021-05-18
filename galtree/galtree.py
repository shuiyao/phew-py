import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import functions
import config
from importlib import reload
import os
from myinit import *

# For each snapshot:
# Each galId maps to a mainId
# Want to know if the mainId is
#  - consumed by any other mainId (root)
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
def load_data():
    model = 'l12n144-phew-movie-200'
    snapnum = 100
    parname = os.path.join(DIRS['DATA'], model, "so_z{:03d}.par".format(snapnum))
    host = pd.read_csv(parname, sep='\s+', names=['galId', 'hostId']).set_index('galId')
    galId2hostId = host.to_dict()['hostId']
    sp = functions.load_stars_from_snapshot(model, snapnum)
    return sp, host

# host:
# galId, hostId: from the .par file
# mainIdHost: The corresponding global mainId of the galaxy with hostId

# gals:


sp, host = load_data()
galId2hostId = host.to_dict()['hostId']
initId2mainId = find_mainId_for_initId(sp)

df=pd.DataFrame({'InitId':initId2mainId.keys(),'mainId':initId2mainId.values()})

gals = sp[['galId','mainId']].drop_duplicates('galId').set_index('galId')
galId2mainId = gals.to_dict()['mainId']
host['mainIdHost'] = host.hostId.map(galId2mainId).astype('Int32')
gals = pd.merge(gals, host.drop('hostId',axis=1), how='left', left_index=True, right_index=True)
gals.mainIdHost.fillna(gals.mainId, inplace=True)




# WANT TO FIND THE DIRECT DESCENDENT OF G0 AT T1

# Galaxy g0: (snapnum=t0, galId=0)
# Galaxy g1: (snapnum=t1>t0, galId=1)

# At time t0, all stars in g0 has the same galId and mainId.
# At time t1, they have different galId(t0) and mainId(t0), but supposedly most of them end up in a single galaxy g0'. 
# If g0.mainId == g0'.mainId, R(g0, g0') = 'SELF'. 
# If g0.mainId <> g0'.mainId, R(g0, g0') = 'MERGE'. 
# 1. Define R(g0, g1) according to the relation between g0' and g1
# g0' can relate to other galId at t1:
# R(g0, g1) = 'SAT' if g0'.hostId = g1
# R(g0, g1) = 'CEN' if g1.hostId = g0'
# R(g0, g1) = 'NGB' if else
# 2. Define R(g0, g1) according to the relation between g0 and g0''
# g0'' at t0 is backtracked from g0':
#   + g0''.mainId = g0'.mainId is found. R(g0, g1) = ''
#   + Not found. R(g0, g1) = 
