#!/bin/bash

model=$1

if [ ! $model ]; then
    echo "Need model name. Force exit."
    exit
fi

pushd select/
bash select.sh $model 1.0
bash select.sh $model 0.2
popd

pushd loadhdf5/
bash loadhdf5.sh $model 1.0
bash loadhdf5.sh $model 0.2
popd

pushd gadget3io/
bash allstars.sh $model
popd

if [ ! -e /scratch/shuiyao/scidata/newwind/$model ]; then
    mkdir /scratch/shuiyao/scidata/newwind/$model    
fi
python match_initwinds_rejoin.py $model z1 078
python match_initwinds_rejoin.py $model z0 100
