#!/bin/bash

HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
DATADIR=$1
snapstrs=(033 058 078 108)
zstrs=(4.0 2.0 1.0 0.0)
snapnums=(33 58 78 108)

while read line
do
    DATADIR=/data002/shuiyao/data/$line
    DESTDIR=$DATADIR/$line/
    # for snapnum in ${snapstrs[@]}; do
    echo mkdir: $DESTDIR
    mkdir $DESTDIR
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
	# 4. Move files to export
	mv $DATADIR/snapshot_${snapstrs[snapi]}.bin $DESTDIR
	mv $DATADIR/snapshot_${snapstrs[snapi]}.aux $DESTDIR
	mv $DATADIR/snapshot_${snapstrs[snapi]}.idnum $DESTDIR
	mv $DATADIR/gal_z${snapstrs[snapi]}.grp $DESTDIR
	mv $DATADIR/gal_z${snapstrs[snapi]}.gtp $DESTDIR
	mv $DATADIR/gal_z${snapstrs[snapi]}.stat $DESTDIR
	mv $DATADIR/so_z${snapstrs[snapi]}.sogrp $DESTDIR
	mv $DATADIR/so_z${snapstrs[snapi]}.sogtp $DESTDIR
	mv $DATADIR/so_z${snapstrs[snapi]}.sovcirc $DESTDIR
	mv $DATADIR/so_z${snapstrs[snapi]}.par $DESTDIR
	# bash skid_one.sh 
    done
#$HDF5EXEC $DATADIR/snapshot_108
done < folders_l12n144.dat
