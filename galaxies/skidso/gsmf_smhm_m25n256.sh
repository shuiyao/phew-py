#!/bin/bash

HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
DATADIR=$1
snapstrs=(033 058 078 108)
zstrs=(4.0 2.0 1.0 0.0)
snapnums=(33 58 78 108)

while read line
do
    DATADIR=/proj/shuiyao/$line
    DESTDIR=$DATADIR/$line/
    # GSMFDIR=/data002/shuiyao/gsmfs_smhms/$line
    GSMFDIR=/scratch/shuiyao/sci/PHEW_TEST/$line
    # for snapnum in ${snapstrs[@]}; do
    echo mkdir: $GSMFDIR
    mkdir $GSMFDIR
    # 4. Generate GSMFs and SMHMs
    echo skid_stats_l25 $line 
    python skid_stats_l25.py $line
    echo smhm_l25 $line 
    python smhm_l25.py $line
#$HDF5EXEC $DATADIR/snapshot_108
done < folders_m25n256.dat
