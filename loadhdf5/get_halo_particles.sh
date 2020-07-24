#!/bin/bash

make
# gdb --args ./get_halo_particles l25n144-phew-fa 78
# ./get_halo_particles l25n144-phew-rcloud 33 25
# ./get_halo_particles l25n288-phew-m5 33 25

# ./get_halo_particles l25n144-phew-rcloud 33 25
# ./get_halo_particles l25n144-phew-rcloud 108 25

# ./get_halo_particles l25n144-phew-m5-spl 108 25
./get_halo_particles l25n144-phew-m5 78 25
./get_halo_particles l25n288-phew-m5 78 25
