#!/bin/tcsh

# 1. Generate .stars file from TIPSY binaries and .grp, .sogrp files
# ./allstars modelname snapnum (all)
# "all" means make no difference on haloes masses

# 2. Generate .starinfo file from the .stars file and SFRINFO files
# python stars_by_mvir modelname snapnum (mh??)

# modelname="p50n288fiducial"
#set modelname="l25n144-phewoff-g3cool"
#set modelname="l25n144-phewoff"
set modelname="l25n144-gadget3"
set flag=0

mkdir /scratch/shuiyao/scidata/gadget3io/$modelname

./allstars $modelname 108
# ./allstars $modelname 58 all
# python stars_by_mvir_gizmo.py $modelname 108 mh11 $flag
# python stars_by_mvir_gizmo.py $modelname 108 mh12 $flag
# python stars_by_mvir_gizmo.py $modelname 108 mh13 $flag
python stars_by_mvir_gadget3.py $modelname 108 mh11 $flag
python stars_by_mvir_gadget3.py $modelname 108 mh12 $flag
python stars_by_mvir_gadget3.py $modelname 108 mh13 $flag
# python allstars.py $modelname 058


