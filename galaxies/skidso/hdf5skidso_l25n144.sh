#!/bin/bash

HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
snapstrs=(033 058 078 108)
zstrs=(4.0 2.0 1.0 0.0)
snapnums=(33 58 78 108)

while read line
do
    DATADIR=/data002/shuiyao/data/$line
    DESTDIR=$DATADIR/$line/
    SCIDIR=/data002/shuiyao/scripts/$line
    GSMFDIR=/data002/shuiyao/gsmfs_smhms/$line
    # for snapnum in ${snapstrs[@]}; do
    echo mkdir: $DESTDIR
    mkdir $DESTDIR
    echo mkdir: $GSMFDIR
    mkdir $GSMFDIR
    echo mkdir: $SCIDIR
    mkdir $SCIDIR
    for snapi in 0 1 2 3; do
        # 1. hdf52tipsy
	echo hdf52tipsy $DATADIR/snapshot_${snapstrs[snapi]}
	$HDF5EXEC $DATADIR/snapshot_${snapstrs[snapi]}
	# 2. skid
	echo skid_one ${snapnums[snapi]} ${zstrs[snapi]} $line
	bash skid_one.sh ${snapnums[snapi]} ${zstrs[snapi]} $line
	# 3. so
	echo so_one ${snapnums[snapi]} ${zstrs[snapi]} $line
	bash so_one.sh ${snapnums[snapi]} ${zstrs[snapi]} $line
    done
    # 4. Generate GSMFs and SMHMs
    echo skid_stats_l25 $line 
    python skid_stats_l25.py $line
    echo smhm_l25 $line
    python smhm_l25.py $line
    echo mzr $line
    bash mzr.sh $line
    echo mrhotz.sh $line
    bash mrhotz.sh $line
#$HDF5EXEC $DATADIR/snapshot_108
done < folders_single.dat
