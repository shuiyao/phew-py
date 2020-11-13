#!/bin/bash

HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
#DATADIR=$1
snapstrs=(033 058 078 108)
zstrs=(4.0 2.0 1.0 0.0)
snapnums=(33 58 78 108)

while read line
do
    echo $line
    DATADIR=/data002/shuiyao/data/$line
    DESTDIR=$DATADIR/$line/
    # REVERT:
    # DATADIR=$DATADIR/$line
    # DESTDIR=/data002/shuiyao/data/$line/

    # for snapnum in ${snapstrs[@]}; do
    echo datadir: $DATADIR
    echo mkdir: $DESTDIR
    mkdir $DESTDIR
    for snapi in 0 1 2 3; do
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
    done
done < folders_single.dat
