#!/usr/bin/bash

#Usage: bash so_one.sh snapnum redshift runname

i=$1
z=$2
runname=$3
SOEXEC=/home/shuiyao/tools/soradl25/sorad
#SOEXEC=/home/shuiyao/tools/soradl12/sorad

database="/proj/shuiyao/"
basename="snapshot_"
gtpbase=$database$runname/"gal_z"
sobase=$database$runname/"so_z"

omega0=0.30
somin=16

if [[ $i -lt 10 ]]
then
    num="00"$i
elif [[ $i -lt 100 ]]
then 
    num="0"$i
else
    num=$i
fi


binfile=$database$runname/$basename$num".bin"
statfile=$gtpbase$num".stat"
sofile=$sobase$num

$SOEXEC $binfile $statfile $sofile $omega0 $z

