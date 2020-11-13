#!/bin/bash

# 1. Generate .stars file from TIPSY binaries and .grp, .sogrp files
# ./allstars modelname snapnum (all)
# "all" means make no difference on haloes masses

# 2. Generate .starinfo file from the .stars file and SFRINFO files
# python stars_by_mvir modelname snapnum (mh??)

modelname=$1
lbox=$2

flag=2
# 0: Gadget3
# 1: GIZMO-PhEW
# 2: GIZMO-PhEW-Extra (Add extra tracking information)

if [ ! $modelname ]; then
    echo "allstars.sh: Need model name. Force exit."
    exit
fi
if [ ! $lbox ]; then
    echo "allstars.sh: Need boxsize. Force exit."
    exit
fi

if [ ! -e /home/shuiyao_umass_edu/scidata/$modelname ]; then
    mkdir /home/shuiyao_umass_edu/scidata/$modelname    
fi

# ./allstars $modelname 78 $lbox all
# python allstars.py $modelname 078 $lbox $flag
# ./allstars $modelname 108 $lbox all
python allstars.py $modelname 108 $lbox $flag

#./allstars $modelname 98 $lbox all
#python allstars.py $modelname 098 $lbox $flag

# Notes:
# ---
# allstars:
#   - Select haloes according to their Msub (from .sogtp)
#   - Returns: All stars selected and their properties
# starinfo_selected_by_....py
#   - load_central_stars() will select only stars from the centrals
#   - find_total_mass() will add up Msub of all stars selected
