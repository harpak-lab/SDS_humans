#!/bin/bash

# Generate job files, ML (total)

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files

# Make file for ML
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/calc_ML_total.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i}.counts -s $OUT/chr${CHR}/chr${CHR}_haps_${i}.snps -w 2 -p 0.1317675 -o $OUT/chr${CHR}/chr${CHR}_haps_${i}" >> do_ml_t_chr${CHR}
done