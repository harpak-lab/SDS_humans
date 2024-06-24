#!/bin/bash

# Generate job files, get results and SEs

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files

# Make file for ML
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/get_boot_pvals_tot.R -m $OUT/chr${CHR}/chr${CHR}_haps_${i}.t.maxl -b $OUT/chr${CHR}/chr${CHR}_haps_${i}.t.boot -o /stor/work/Kirkpatrick/scratch/Jared/ukb/analysis_new_hap/hap_data/chr${CHR}/chr${CHR}_haps_${i}" >> do_getresult_tot_chr${CHR}
done