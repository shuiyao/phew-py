#!/bin/bash

fbase="/proj/shuiyao/"
simname="l25n144phew"
lbox="25"
# fbase="/scratch/shuiyao/data/"
# simname="p25n144rwlIII"
# lbox="25"


cd $fbase$simname/SFRINFO/
if ! [[ -e c.sfrinfo ]]
then
    echo "Create c.sfrinfo ... "
    cat sfrinfo.* > c.sfrinfo
    ls -lh c.sfrinfo
else
    echo "c.sfrinfo exists."
    ls -lh c.sfrinfo    
fi
cd -
# accsingle.py: Generation accretion info
# Input: $SFRINFO/c.sfrinfo
# Output: $SCIDATA/newwind/acc.$simname
#   - z, cold, hot, wind, windret, zcold, zhot, zwind
python accsingle.py $lbox $fbase $simname 1

# flag = 1: PhEW
# flag = 0: Non-PhEW

# plotacc.py: Plot the accretion info from multiple simulations (def. 2)
