#!/bin/bash

# mode = 0: sort .sogrp data and show
# mode = 1: Run get_halo_particles.c

mode=$1

# datadir=/proj/shuiyao/l25n144-phewoff
datadir=/proj/shuiyao/p50n576fi
# snapname=$datadir/snapshot_108
# sogrpname=$datadir/so_z108.sogrp
# sovcircname=$datadir/so_z108.sovcirc
snapname=$datadir/snap_p50n576zw_108
sogrpname=$datadir/gals/so_z108.sogrp
sovcircname=$datadir/gals/so_z108.sovcirc

if [ $mode -eq 0 ]; then
    cp $sovcircname ./temp.dat
    sort -g -r -k2,2 temp.dat > sostat.dat
    rm -f ./temp.dat
fi

if [ $mode -eq 1 ]; then
    ./get_halo_particles $snapname $sogrpname 13166
    ./get_halo_particles $snapname $sogrpname 17125
    ./get_halo_particles $snapname $sogrpname 18476
    ./get_halo_particles $snapname $sogrpname 9554
    ./get_halo_particles $snapname $sogrpname 13770
    ./get_halo_particles $snapname $sogrpname 14228
    ./get_halo_particles $snapname $sogrpname 19580            
fi
