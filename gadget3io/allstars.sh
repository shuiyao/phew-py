#!/bin/bash

# 1. Generate .stars file from TIPSY binaries and .grp, .sogrp files
# ./allstars modelname snapnum (all)
# "all" means make no difference on haloes masses

# 2. Generate .starinfo file from the .stars file and SFRINFO files
# python stars_by_mvir modelname snapnum (mh??)

modelname=$1
lbox=$2
ncpu=$3

flag=2
# 0: Gadget3
# 1: GIZMO-PhEWOff
# 2: GIZMO-PhEW-Extra (Add extra tracking information)

if [ ! $modelname ]; then
    echo "allstars.sh: Need model name. Force exit."
    exit
fi
if [ ! $lbox ]; then
    echo "allstars.sh: Need boxsize. Force exit."
    exit
fi
if [ ! $ncpu ]; then
    ncpu=256
    echo "Use default ncpu = 256"
fi

if [ ! -e /home/shuiyao_umass_edu/scidata/$modelname ]; then
    mkdir /home/shuiyao_umass_edu/scidata/$modelname    
fi

#./allstars $modelname 108 $lbox all
# .stars
#Idx ID GID HID Mass Tmax Age
#46021453    1633394     0  5265 1.17962e-09 4.978 0.66032

python allstars.py $modelname 108 $lbox $ncpu $flag
# .starinfo:
#a_form a_acc a_last Mass WindMass StarMass Tmax Z GID HID 
#0.42882 0.37561  0.00000 1.85705e+41 1.14229e+41 0.00000e+00  4.849 0.0 2610 2610

# Notes:
# ---
# allstars:
#   - Select haloes according to their Msub (from .sogtp)
#   - Returns: All stars selected and their properties
# starinfo_selected_by_....py
#   - load_central_stars() will select only stars from the centrals
#   - find_total_mass() will add up Msub of all stars selected
