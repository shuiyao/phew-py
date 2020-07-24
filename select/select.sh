#!/bin/bash

# simname="m6n64beta6"
# ncpu="32"
# simname="l25n144-phew"
# simname="l25n144-phew-m4kh100fs10"
simname=$1
redshift=$2
ncpu="256"

if [ ! $simname ]; then
    echo "select.sh: Need model name. Force exit."
    exit
fi
if [ ! $redshift ]; then
    echo "Use default redshift = 1.0"
    redshift="1.0"
fi
if [[ $redshift == "2.0" ]]; then
    zstr="z2/"
elif [[ $redshift == "1.0" ]]; then
    zstr="z1/"
elif [[ $redshift == "0.2" ]]; then
    zstr="z0/"
fi

fbase="/proj/shuiyao/"$simname"/WINDS/"
odir=$fbase$zstr
mkdir $odir
# Input: $WINDS/winds.*
# Output: $WINDS/z?/wid.dat; an ID array for all wind particles launched within the time window

echo "---- Selecting WIND files ----"
python select_winds.py $ncpu $redshift $fbase

# Input: $WINDS/phews.* (The format for PhEW and GW are different); $WINDS/z?wid.dat
# Output: $WINDS/z?/phews.*
echo "---- Selecting REJOIN files ----"
python select_rejoin.py $ncpu $redshift $fbase
echo "---- Selecting PHEWTRACK files ----"
python select_phewtracks.py $ncpu $redshift $fbase

# echo "---- Selecting TRACK files ----"
#python select_tracks.py $ncpu $redshift $fbase

python find_host_haloes_for_phews.py $ncpu $redshift $simname

cd $odir
# Prepare initwinds.sorted and rejoin.sorted
rm -f initwinds.sorted
cat initwinds.* > initwinds.all
sort -n -k2,2 initwinds.all > initwinds.sorted
rm -f initwinds.all

rm -f rejoin.sorted
cat rejoin.* > rejoin.all
sort -n -k2,2 rejoin.all > rejoin.sorted
rm -f rejoin.all
# cat track.* > track.all
# sort -n -k2,2 track.all > sorted.tracks
# rm -f track.all

# Prepare phews.sorted
# SHOULD CHECK HERE IF THE FILES IS TOO BIG

# sizetot=`du | awk -F " " '{print $1}'`
# echo "Total Size = "$sizetot" KB"
# if [[ $sizetot -gt "524288" ]] # 512M
# then
#     echo "Total Size > 524288, Exit!"
#     exit 110
# fi

du -h
rm -f sorted.phews
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
