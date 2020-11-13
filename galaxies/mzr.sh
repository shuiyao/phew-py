#!/bin/bash

CODE=/scratch/shuiyao/sci/mzr/gal_metals
DATA_FOLDER=/proj/shuiyao
#MODEL=l25n144-phewoff
MODEL=l25n144-gizmo-mufasa
#MODEL=l25n144-g3coolsf

$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 108 $MODEL/mzr_108.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 78 $MODEL/mzr_078.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 58 $MODEL/mzr_058.txt
$CODE $DATA_FOLDER/$MODEL/snapshot_ $DATA_FOLDER/$MODEL/gal_z 33 $MODEL/mzr_033.txt
