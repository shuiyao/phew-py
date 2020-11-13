#!/bin/bash

# simname="p25n144rwlII"
# fbase="/scratch/shuiyao/data/"$simname"/WINDS/"
simname="p25n144rwlIII"
fbase="/scratch/shuiyao/data/"$simname"/WINDS/"
# simname="p25n144phewV"
# fbase="/work/shuiyao/data/"$simname"/WINDS/"
redshift="0.2"
ncpu="128"

python select_winds_by_mvir.py $ncpu $redshift $fbase
python select_tracks_by_mvir.py $ncpu $redshift $fbase
