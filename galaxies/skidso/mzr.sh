#!/bin/bash

CODE=/scratch/shuiyao/sci/mzr/gal_metals
DATA_FOLDER=/proj/shuiyao
MODEL=$1
OUT_FOLDER=/scratch/shuiyao/sci/PHEW_TEST/$MODEL

# $CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 108 $MODEL/mzr_108.txt
# $CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 78 $MODEL/mzr_078.txt
# $CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 58 $MODEL/mzr_058.txt
# $CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 33 $MODEL/mzr_033.txt

$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 108 $OUT_FOLDER/mzr_108.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 78 $OUT_FOLDER/mzr_078.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 58 $OUT_FOLDER/mzr_058.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 33 $OUT_FOLDER/mzr_033.txt

