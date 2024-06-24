#!/bin/bash

# Generate job files, get results and SEs

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files

# Make file for ML
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/get_boot_pvals_fec.R -m $OUT/chr${CHR}/chr${CHR}_haps_${i}.f.maxl -b $OUT/chr${CHR}/chr${CHR}_haps_${i}.f.boot -o $OUT/chr${CHR}/chr${CHR}_haps_${i}" >> do_getresult_fec_chr${CHR}
done