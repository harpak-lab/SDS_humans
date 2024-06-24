#!/bin/bash

# Generate job files, transpose

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files
let ENDNUM=NUM-1

# Make file for ML

#chromosome chunks

for i in $(seq 0 $NUM); do
        echo "Rscript /scratch/ukb/scripts/Transpose_haplotypes.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i} -s $OUT/chr${CHR}/chr${CHR}_haps_lrs.sample -v $OUT/chr${CHR}/chr${CHR}.snps.txt -o $OUT/chr${CHR}/chr${CHR}_haps_${i}" >> do_transpose_chr${CHR}
done