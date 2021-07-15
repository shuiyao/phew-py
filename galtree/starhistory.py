import matplotlib.pyplot as plt
import pandas as pd
import functions
import config
from importlib import reload
import os

from myinit import *

reload(functions)

# model = 'l25n144-phew-m5-spl'
model = 'l12n144-phew-movie-200'

def generate_star_history(model, start=0, skip=1):
    '''
    Trace the location of each star over cosmic time and assign the mainIds for each star at different snapshots. One can generate the merger history of galaxies from the outputs.
    Will write a CSV file for each snapshot explored. The CSV file contains all star particles within that snapshot, with the following properties:
    sid:    The unique particle Id that is constant over snapshots.
    galId:  The ID of the host galaxy at that snapshot.
    mainId: A unique identifier for galaxies that is constant across snapshots.
    initId: The mainId when the star first appeared in a galaxy, constant or 0.
    mass:   The mass of the star at the snapshot.

    Parameters
    ----------
    model: str.
        The name of the model
    start: int, default 0.
        The starting snapshot to do the tracking.
        If non-zero, restart from existing output at that snapshot.
    '''

    path_model = os.path.join(DIRS['DATA'], model)
    if(start == 0):
        config.spAll = pd.DataFrame(columns=['sid','galId','mainId','initId','mass']).set_index('sid')
        config.maxMainId = 0
        for i in range(201)[::skip]:
            functions.process_snapshot(model, i)
    else:
        # Restart from some snapshots
        config.spAll = functions.load_stars_from_snapshot(path_model, start)
        config.maxMainId = config.spAll.mainId.max()
        for i in range(201)[start+skip::skip]:
            functions.process_snapshot(model, i)
    



    
