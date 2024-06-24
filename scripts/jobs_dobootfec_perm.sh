#!/bin/bash

# Generate job files, bootstrapping (permuted)
# fecundity

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files


# Make file for bootstrapping
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/calc_boot_fecundity.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i}.permuted.counts -w 2 -p 0.1317675 -b 100 -m $OUT/chr${CHR}/chr${CHR}_haps_${i}.malematrix.permuted.rds -g $OUT/chr${CHR}/chr${CHR}_haps_${i}.femalematrix.permuted.rds -o $OUT/chr${CHR}/chr${CHR}_haps_${i}.permuted" >> do_permboot_fec_chr${CHR}
done