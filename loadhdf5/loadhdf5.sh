#!/bin/bash

folder="/proj/shuiyao/"
# modelname="l25n144-phew-m4kh100fs10"
modelname=$1
redshift=$2

if [ ! $modelname ]; then
    echo "loadhdf5.sh: Need model name. Force exit."
    exit
fi
if [ ! $redshift ]; then
    echo "Use default redshift = 1.0"
    redshift="1.0"
fi
if [[ $redshift == "1.0" ]]; then
    snapstr="078"
    snapnum=78
elif [[ $redshift == "0.25" ]]; then
    snapstr="098"
    snapnum=98
elif [[ $redshift == "0.2" ]]; then
    snapstr="100"
    snapnum=100
fi

fbase=$folder$modelname"/"
snap=$fbase"snapshot_"$snapstr

echo $snap

./get_particles -phew -sph $snap
./rhot $fbase"snapshot" z0 $snapnum 40 40
mv mrhot_phew_z0_$snapstr /scratch/shuiyao/sci/PHEW_TEST/$modelname/
./rhot $fbase"snapshot" z0 $snapnum
mv mrhot_phew_z0_$snapstr /scratch/shuiyao/sci/PHEW_TEST/$modelname/mrhot_z0_$snapstr
# gdb --args ./loadhdf5 $snap
