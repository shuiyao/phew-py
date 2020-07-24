#!/bin/bash

# 1. Generate .stars file from TIPSY binaries and .grp, .sogrp files
# ./allstars modelname snapnum (all)
# "all" means make no difference on haloes masses

# 2. Generate .starinfo file from the .stars file and SFRINFO files
# python stars_by_mvir modelname snapnum (mh??)

# modelname="p50n288fiducial"
#set modelname="l25n144-phewoff-g3cool"
#set modelname="l25n144-phewoff"
#set modelname="l25n144-gadget3"
#set modelname="l25n144-phewoff"

modelname=$1
lbox=$2
flag=1 # PHEW FLAG: 1 = PHEW

if [ ! $modelname ]; then
    echo "allstars.sh: Need model name. Force exit."
    exit
fi
if [ ! $lbox ]; then
    echo "allstars.sh: Need boxsize. Force exit."
    exit
fi

# !!!!! DON'T FORGET TO CHANGE UNIT_M in allstars.py
mkdir /scratch/shuiyao/scidata/gadget3io/$modelname

./allstars $modelname 88 $lbox all
python allstars.py $modelname 088 $lbox $flag

# ./allstars $modelname 108
# python stars_by_mvir_gizmo.py $modelname 108 mh11 $flag
# python stars_by_mvir_gizmo.py $modelname 108 mh12 $flag
# python stars_by_mvir_gizmo.py $modelname 108 mh13 $flag
# python stars_by_mvir_gadget3.py $modelname 108 mh11 $flag
# python stars_by_mvir_gadget3.py $modelname 108 mh12 $flag
# python stars_by_mvir_gadget3.py $modelname 108 mh13 $flag
# python allstars.py $modelname 058


# Notes:
# ---
# allstars:
#   - Select haloes according to their Msub (from .sogtp)
#   - Returns: All stars selected and their properties
# stars_by_mvir_gizmo.py
#   - Combine output above with SFRINFO
#   - 
# starinfo_selected_by_....py
#   - load_central_stars() will select only stars from the centrals
#   - find_total_mass() will add up Msub of all stars selected
