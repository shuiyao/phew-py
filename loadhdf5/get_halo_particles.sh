#!/bin/bash

model=$1

if [ ! $model ]; then
    echo "get_halo_particles.sh: Need model name."
    exit
fi    

#make
# gdb --args ./get_halo_particles l25n144-phew-fa 78
# ./get_halo_particles l50n288-phew-m5 98 50
# ./get_halo_particles l25n288-phew-m5 98 25
# ./get_halo_particles l50n288-phew-m4 98 50
# ./get_halo_particles l50n288-phew-m5 78 50
# ./get_halo_particles l25n288-phew-m5 78 25
# ./get_halo_particles l50n288-phew-m4 78 50


# Turn off PHEW_EXTRA_OUTPUT
./get_halo_particles $model 98 50
