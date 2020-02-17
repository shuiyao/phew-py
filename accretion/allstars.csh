#!/bin/tcsh

# 1. Generate .stars file from TIPSY binaries and .grp, .sogrp files
# ./allstars modelname snapnum (all)
# "all" means make no difference on haloes masses

# 2. Generate .starinfo file from the .stars file and SFRINFO files
# python stars_by_mvir modelname snapnum (mh??)

modelname="p50n288fiducial"

./allstars $modelname 108
./allstars $modelname 58 all
python stars_by_mvir.py $modelname 108 mh11
python stars_by_mvir.py $modelname 108 mh12
python stars_by_mvir.py $modelname 108 mh13
python allstars.py $modelname 058


