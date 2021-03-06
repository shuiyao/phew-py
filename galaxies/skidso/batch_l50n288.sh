#!/bin/bash

#HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
HDF5EXEC=/home/shuiyao/tools/SimIO/hdf52tipsy/hdf52tipsy
NOSKID=$1
snapstrs=(033 058 078 108 098)
zstrs=(4.0 2.0 1.0 0.0 0.25)
snapnums=(33 58 78 108 98)

while read line
do
    DATADIR=/proj/shuiyao/$line
    DESTDIR=$DATADIR/$line/
    # for snapnum in ${snapstrs[@]}; do
    echo mkdir: $DESTDIR
    mkdir $DESTDIR
    mkdir /scratch/shuiyao/sci/PHEW_TEST/$line
    if [ ! $NOSKID ]; then
	for snapi in 0 1 2 3 4; do
            # 1. hdf52tipsy
	    echo hdf52tipsy $DATADIR/snapshot_${snapstrs[snapi]}
	    $HDF5EXEC $DATADIR/snapshot_${snapstrs[snapi]}
	    # 2. skid
	    echo skid_one ${snapnums[snapi]} ${zstrs[snapi]} $line
	    bash skid_one.sh ${snapnums[snapi]} ${zstrs[snapi]} $line
	    # 3. so
	    echo so_one ${snapnums[snapi]} ${zstrs[snapi]} $line
	    bash so_one_l50n288.sh ${snapnums[snapi]} ${zstrs[snapi]} $line
	done
    fi
    # 4. GSMFs and SMHMs
    # echo gsmf_smhm.sh $line
    # bash gsmf_smhm.sh $line
    # 5. Metal Rho-T Diagram
    # echo mrhotz.sh $line
    # bash mrhotz.sh $line
    # 6. MZR and fgas
    # echo mzr.sh $line
    # bash mzr.sh $line
    echo fgas.sh $line
    bash fgas.sh $line
#$HDF5EXEC $DATADIR/snapshot_108
done < folders_l50n288.dat
