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
elif [[ $redshift == "0.2" ]]; then
    snapstr="100"
fi

fbase=$folder$modelname"/"
snap=$fbase"snapshot_"$snapstr

echo $snap

./loadhdf5 -phew -sph $snap
# gdb --args ./loadhdf5 $snap
