#!/bin/bash

HDF5EXEC=/home/shuiyao_umass_edu/SimIO/hdf52tipsy/hdf52tipsy
model=$1
lbox=$2
NOSKID=$3
snapstrs=(033 058 078 108 098 088)
zstrs=(4.0 2.0 1.0 0.0 0.25 0.5)
snapnums=(33 58 78 108 98 88)

DATADIR=/nas/astro-th-nas/shuiyao/$model
DESTDIR=/home/shuiyao_umass_edu/scidata/$model
mkdir $DESTDIR
if [ ! $NOSKID ]; then
    for snapi in 0 1 2 3 4; do
        # 1. hdf52tipsy
	echo hdf52tipsy $DATADIR/snapshot_${snapstrs[snapi]}
	$HDF5EXEC $DATADIR/snapshot_${snapstrs[snapi]}
	# 2. skid
	echo skid_one ${snapnums[snapi]} ${zstrs[snapi]} $model
	bash skid_one.sh ${snapnums[snapi]} ${zstrs[snapi]} $model
	# 3. so
	echo so_one ${snapnums[snapi]} ${zstrs[snapi]} $model
	bash so_one_l25n288.sh ${snapnums[snapi]} ${zstrs[snapi]} $model
    done
fi
# 4. Generate GSMFs and SMHMs
echo "   - Compiling Galaxies ..."
python ./scripts/gsmfs.py $model $lbox
echo "   - Compiling Haloes ..."
python ./scripts/smhms.py $model $lbox

# 5. Metal Rho-T Diagram
# echo "Calculating Phase Diagrams"
# for snapi in 0 1 2 3 4; do
#     ../loadhdf5/rhot $DATADIR/snapshot ${snapnums[snapi]}
#     mv mrhot_phew_${snapstrs[snapi]} $DESTDIR/
# done

# 6. MZR and fgas
# echo "Calculating Mass-Metallicity Relation..."
# bash mzr.sh $model
# echo "Calculating Gas Fractions..."
# bash fgas.sh $model

