#!/bin/bash

# Generate job files, ML permuted (fecundity)

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files


# Make file for ML
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/calc_ML_fecundity.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i}.permuted.counts -s $OUT/chr${CHR}/chr${CHR}_haps_${i}.snps -w 2 -p 0.1317675 -o $OUT/chr${CHR}/chr${CHR}_haps_${i}.permuted" >> do_ml_f_perm_chr${CHR}
done