#!/bin/bash

folder="/nas/astro-th-nas/shuiyao/"
modelname=$1
redshift=$2

if [ ! $modelname ]; then
    echo "loadhdf5.sh: Need model name. Force exit."
    exit
fi
if [ ! $redshift ]; then
    echo "Use default redshift = 0.25"
    redshift="0.25"
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

./get_particles -phew $snap

# Output:
# $datadir/$model/snapshot_???.phews
# idx rhoc LogTc rhoa LogTa fc fw LastSFTime

./rhot $fbase"snapshot" $snapnum
mv mrhot_phew_$snapstr /home/shuiyao_umass_edu/sci/phew-py/data/$modelname/

# Output:
# ncells_x, ncells_y
# xnodes
# ynodes
# number_of_particles total_mass wind_mass

# The default grid size is 256 x 256
# Although we can have a smaller grid, e.g., 40 x 40
# ./rhot $fbase"snapshot" z0 $snapnum 40 40


# gdb --args ./loadhdf5 $snap
