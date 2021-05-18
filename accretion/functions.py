import numpy as np
import pandas as pd

def load_timeinfo_for_snapshots():
    df_redz = pd.read_csv("../redshifts.txt", sep='\s+', header=0)
    df_redz = df_redz.drop("#snapnum", axis=1)
    return df_redz

# Find the first snapshot after ainit
snapnum_first = (df_redz.a >= ainit).argmax()

# Select all PhEW particles inside a given halo at snapnum
gpTable.query("snapnum==@snapnum AND galId==galId AND Mc > 0")

def gasp_mass_gain_since_last_snapshot(PID, snapshot, gpTable):
    mass_this = gpTable.query('snapnum==@snapnum AND PId==@PID').Mwind
    mass_last = gpTable.query('snapnum==@snapnum-1 AND PId==@PID').Mwind
    return mass_this - mass_last

def phew_mass_loss_since_last_checkpoint(PID, snapnum, gpTable, initTable):
    '''
    Find the amount of mass loss for a PhEW particle (PID) between the 
    last checkpoint MAX(a{snapnum-1}, ainit)

    Parameters
    ----------
    PID: int
         The Particle ID of the PhEW particle
    snapnum: int
         The end time for query.
    gpTable: 
         The table for all gas particles
    initTable:
         The table for the initwinds information

    Return
    ------
    Mass loss of the PhEW particle during the given period of time.
    '''

    # gpTable.type == pd.DataFrame:
    mass_this = gpTable.query('snapnum==@snapnum AND PId==@PID').Mass
    ainit = phewTable.query('PId==@PID').ainit
    minit = phewTable.query('PId==@PID').minit
    if(redzTable.loc[snapnum-1].a < ainit):
        mass_last = minit
    else:
        mass_last = gpTable.query('snapnum==@snapnum-1 AND PId==@PID').Mass
    return mass_this - mass_last

# Mass Change From Last Checkpoint
# 1. Normal Gas Particle: dMass = Mwind[i] - Mwind[i-1] or Mwind[acc] - Mwind[i-1]
#    (snapnum, PID) -> Mgain
# 2. PhEW Particle: dMass = Mass[i] - Mass[i-1] or Mass[i] - Mass[init]
#    (snapnum, PID) -> Mloss
