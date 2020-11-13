#!/bin/bash

CODE=/scratch/shuiyao/sci/fgas/fgas_mstar
DATA_FOLDER=/proj/shuiyao
MODEL=$1
OUT_FOLDER=/scratch/shuiyao/sci/PHEW_TEST/$MODEL

$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/ 108 $OUT_FOLDER/fgas_mstar_108.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/ 78 $OUT_FOLDER/fgas_mstar_078.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/ 58 $OUT_FOLDER/fgas_mstar_058.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/ 33 $OUT_FOLDER/fgas_mstar_033.txt

