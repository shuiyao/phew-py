import matplotlib.pyplot as plt
import pandas as pd
import functions
import config
from importlib import reload
import os

reload(functions)

# model = 'l25n144-phew-m5-spl'
model = 'l12n144-phew-movie-200'

def generate_star_history(model, start=0):
    if(start == 0):
        config.spAll = pd.DataFrame(columns=['sid','galId','mainId','initId','mass']).set_index('sid')
        config.maxMainId = 0
        for i in range(200)[::5]:
            functions.process_snapshot(model, i)
    else:
        # Restart from some snapshots
        config.spAll = functions.load_stars_from_snapshot(model, start)
        config.maxMainId = config.spAll.mainId.max()
        for i in range(200)[start+5::5]:
            functions.process_snapshot(model, i)



    
