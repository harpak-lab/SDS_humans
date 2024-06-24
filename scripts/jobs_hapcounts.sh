#!/bin/bash

# Create haplotype counts submission for local cluster

OUT=/scratch/ukb/data/	

#Arguments
CHR=$1
NUM=$2  #number of files

#make file
for i in $(seq 0 $NUM); do
	echo "Rscript /scratch/ukb/scripts/haplotype_counts.R -f $OUT/chr${CHR}/chr${CHR}_haps_${i}.hap -w 2 -c ${CHR} -v $OUT/chr${CHR}/chr${CHR}_haps_${i}.snps -o $OUT/chr${CHR}/chr${CHR}_haps_${i}" >> do_hapcounts_chr${CHR}
done
