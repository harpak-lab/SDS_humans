#!/bin/bash

# Generate job files, bootstrapping (observed)

OUT=/scratch/ukb/data/

#Arguments
CHR=$1
NUM=$2  #number of files

# Make file for bootstrapping
#chromosome chunks

for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/calc_boot_viability.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i}.counts -w 2 -p 0.1317675 -b 100 -o $OUT/chr${CHR}/chr${CHR}_haps_${i}" >> do_boot_chr${CHR}
	 done