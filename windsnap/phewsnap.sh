#!/bin/bash

fbase="/proj/shuiyao/"

modelname="l25n144-phew-rcloud"
#modelname="l25n144-phewoff"
binfile=$fbase$modelname"/snapshot_098.bin"
auxfile=$fbase$modelname"/snapshot_098.aux"
awfile=$fbase$modelname"/snapshot_098.aw"

mkdir /scratch/shuiyao/scidata/windsnap/$modelname

echo $binfile
echo $auxfile
echo $modelname
echo $awfile
./phewsnap $binfile $auxfile $modelname $awfile -0.1 0.1


