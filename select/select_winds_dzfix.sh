#!/bin/bash

# simname="m6n64beta6"
# ncpu="32"
# simname="l25n144-phew"
# simname="l25n144-phew-m4kh100fs10"
simname=$1
ncpu=256

fbase="/nas/astro-th/shuiyao/"$simname"/WINDS/"

python select_winds_dzfix.py $ncpu 4.0 $fbase
python select_winds_dzfix.py $ncpu 2.0 $fbase
python select_winds_dzfix.py $ncpu 1.0 $fbase
python select_winds_dzfix.py $ncpu 0.2 $fbase
