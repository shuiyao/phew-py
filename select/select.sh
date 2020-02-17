#!/bin/bash

simname="m6n64beta6"
ncpu="32"
# simname="l25n144phew"
# ncpu="128"

fbase="/proj/shuiyao/"$simname"/WINDS/"
odir=$fbase"z2/"
redshift="2.0"
mkdir $odir
# Input: $WINDS/winds.*
# Output: $WINDS/z?/wid.dat; an ID array for all wind particles launched within the time window
python select_winds.py $ncpu $redshift $fbase
# Input: $WINDS/phews.* (The format for PhEW and GW are different); $WINDS/z?wid.dat
# Output: $WINDS/z?/phews.*
python select_tracks.py $ncpu $redshift $fbase

cd $odir
# SHOULD CHECK HERE IF THE FILES IS TOO BIG
sizetot=`du | awk -F " " '{print $1}'`
echo "Total Size = "$sizetot" KB"
if [[ $sizetot -gt "524288" ]] # 512M
then
    echo "Total Size > 524288, Exit!"
    exit 110
fi
du -h
cat phews.* > all.phews
sort -n -k2,2 all.phews > sorted.phews
rm -f all.phews

cd -
# Most relevant for PhEW:
# select_winds_by_mvir.py: Select winds from certain time window, and group them according to their host halo mass
# Input: $SIM/pid_z*.grp, $WINDS/winds.*
# Output: $WINDS/z*/wid???.dat
#  - wid105: 10.3 < mvir < 10.7
#  - wid110: 10.9 < mvir < 11.1
#  - wid115: 11.4 < mvir < 11.6
#  - wid120: 11.9 < mvir < 12.1
#  - wid125: 12.5 < mvir < 13.0
# python select_winds_by_mvir.py $ncpu $redshift $fbase

# Input: $WINDS/z*/sorted.tracks, $WINDS/z*/wid???.dat
# Output: $WINDS/z*/tracks???
# python select_tracks_by_mvir.py $ncpu $redshift $fbase
