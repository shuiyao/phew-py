#!/bin/bash

# simname="p50n288fofVIII"
while read simname; do
    fbase=/data002/shuiyao/data/$simname/WINDS/
    ncpu=256
    python pickout_winds_4z.py $ncpu $fbase
    echo $fbase
#done < p25n144XIII.lst
#done < p50n288fof.lst
done < temp.lst

