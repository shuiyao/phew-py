# Working directly with the sfrinfo files

from myinit import *
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

model = "l50n288-phew-m5"
#model = "l25n288-phew-m5"
#model = "l50n288-phewoff-old"
#model = "p50n288fiducial"

fname = DIRS['DATA'] + model + "/SFRINFO/sfrinfo.all"

tab = pd.read_csv(fname, sep='\s+', header=0)
tab.rename(columns=({'#a':'a'}), inplace=True)

sub = tab[tab['LastSFTime'] > 0]
normal = tab[tab['LastSFTime'] == 0]
subw = tab[tab['LastSFTime'] < 0]

print ("Number of spurious accretion: %d/%d" % (sub.shape[0], tab.shape[0]))
# da = sub['a'] - sub['LastSFTime']
# from cosmology import tcosmic
# dt = (tcosmic(sub['a']) - tcosmic(sub['LastSFTime'])) / 1.e6
# dt95 = np.percentile(dt, 95)
# dt50 = np.percentile(dt, 50)
