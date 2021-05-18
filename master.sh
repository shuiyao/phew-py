#!/bin/bash

modelname=$1

if [ ! $modelname ]; then
    echo "master.sh: Need to specify modelname."
    echo "Usage: bash master.sh modelname"    
fi

flog=$modelname.masterlog

pushd loadhdf5/

echo "loadhdf5." >> flog
bash loadhdf5.sh $modelname 0.25
# - get_particles -phew $snap
# - rhot $fbase"snapshot" z0 $snapnum 256 256

echo "get_halo_particles." >> flog
./get_halo_particles $modelname 98 50

echo "calc_halogas_components." >> flog
python calc_halogas_components $modelname 58
python calc_halogas_components $modelname 108
# Input: gal, so, hdf5
# default lbox = 50, Tcut = 1.e5 K

echo "calc_radial_profile." >> flog
#calc_radial_profile.py

popd

pushd gadget3io/
bash allstars.sh $modelname 50 1024
#bash allstars_phewoff.sh

# show_halogas_components.py
# show_radial_profile.py
# show_phase_diagrams.py
