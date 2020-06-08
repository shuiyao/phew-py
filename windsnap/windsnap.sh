#!/bin/bash

fbase="/proj/shuiyao/"

# modelname="l25n144-phew-rcloud"
modelname="l25n144-phewoff"
binfile=$fbase$modelname"/snapshot_098.bin"
auxfile=$fbase$modelname"/snapshot_098.aux"

mkdir /scratch/shuiyao/scidata/windsnap/$modelname

echo $binfile
echo $auxfile
echo $modelname

./windsnap $binfile $auxfile $modelname -0.1 0.1


