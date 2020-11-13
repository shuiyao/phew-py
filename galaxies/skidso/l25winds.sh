#!/bin/bash

# Do select.sh first
# Purpose: get windsinfo files, ID Mvir Rvir Vinit V25 Rret

fbase="/data002/shuiyao/data/"
redshift="2.0" # <--------------------------------
#simname="l25n144-phewoff" # <--------------------------------
simname="l25n144-gadget3" # <--------------------------------
mkdir windsinfo/$simname
cd $fbase$simname
idnum2ascii snapshot_058 # <--------------------------------
cd - 
python l25winds.py $redshift $fbase $simname # <--------------------------------


