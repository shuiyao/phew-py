#!/bin/bash

#HDF5EXEC=/home/shuiyao/SimIO/hdf52tipsy/hdf52tipsy
HDF5EXEC=/home/shuiyao/tools/SimIO/hdf52tipsy/hdf52tipsy
DATADIR=$1
snapstrs=(033 058 078 108 100)
zstrs=(4.0 2.0 1.0 0.0 0.2)
snapnums=(33 58 78 108 100)

while read line
do
    DATADIR=/proj/shuiyao/$line
    DESTDIR=$DATADIR/$line/
    # for snapnum in ${snapstrs[@]}; do
    i=0
    while read zline
    do
	if [[ $i -ge 0 && $i -lt 109 ]]; then
	    num=$i
	    if [ ${#num} -eq 1 ]
	    then
		num="00"$num
	    elif [ ${#num} -eq 2 ]
	    then
		num="0"$num
	    elif [ ${#num} -eq 3 ]
	    then
		num=$num
	    else
		echo "Wrong Parameter."
		exit 0
	    fi
            # 1. hdf52tipsy
	    echo hdf52tipsy $DATADIR/snapshot_$num
	    $HDF5EXEC $DATADIR/snapshot_$num
	    # 2. skid
	    echo skid_one $i $zline $line
	    bash skid_one.sh $i $zline $line
	    # 3. so
	    echo so_one $i $zline $line
	    bash so_one.sh $i $zline $line
	fi
	i=$(($i+1))	
    done < "redshifts_p50n288ezw15.txt"
done < folders_skidall.dat
