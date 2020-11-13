#!/bin/bash

# Do l25winds.sh afterwards

simname="l25n144-gadget3"
fbase="/data002/shuiyao/data/"$simname"/WINDS/"
odir=$fbase"z2/"
redshift="2.0"
ncpu="128"
mkdir $odir
python select_winds_l25.py $ncpu $redshift $fbase
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
cat track.* > all.tracks
sort -n -k2,2 all.tracks > sorted.tracks
rm -f all.tracks

cd -
#python select_winds_by_mvir.py $ncpu $redshift $fbase
#python select_tracks_by_mvir.py $ncpu $redshift $fbase
