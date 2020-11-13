#!/bin/bash -x

CODE_MRHOT=/scratch/shuiyao/sci/rhot/rhot_mhist
CODE_METAL=/scratch/shuiyao/sci/metals/metals
DATA_FOLDER=/proj/shuiyao
#MODEL=l25n144-phew-mach1
#MODEL=l25n144-phewoff
#MODEL=l25n144-gadget3
MODEL=l25n144-phew-fa

# echo $CODE_MRHOT $DATA_FOLDER/$MODEL/snapshot z0 100
# $CODE_MRHOT $DATA_FOLDER/$MODEL/snapshot z0 100
# mv mrhot_z0_100 ./$MODEL/

echo $CODE_METAL $DATA_FOLDER/$MODEL/snapshot z0 98
$CODE_METAL $DATA_FOLDER/$MODEL/snapshot z0 98
mv tabion_z0_098 ./$MODEL/

# echo $CODE_MRHOT $DATA_FOLDER/$MODEL/snapshot z1 78
# $CODE_MRHOT $DATA_FOLDER/$MODEL/snapshot z1 78
# echo $CODE_METAL $DATA_FOLDER/$MODEL/snapshot z1 78
# $CODE_METAL $DATA_FOLDER/$MODEL/snapshot z1 78
# mv mrhot_z1_078 ./$MODEL/
# mv tabmet_z1_078 ./$MODEL/

