#!/bin/bash

# Show what kinds of things have been done for a particular model

model=$1

SCI=/scratch/shuiyao/sci
SCIDATA=/scratch/shuiyao/scidata
DATA=/proj/shuiyao

echo $DATA/$model/WINDS/":"
dir $DATA/$model/WINDS/

echo $SCI/PHEW_TEST/$model":"
ls $SCI/PHEW_TEST/$model

echo ""
echo $SCIDATA/newwind/$model":"
ls $SCIDATA/newwind/$model

echo ""
echo /proj/shuiyao/$model/*.phews":"
ls /proj/shuiyao/$model/*.phews

