#!/bin/bash

# Create haplotype counts submission for local cluster

#Arguments
CHR=$1
NUM=$2  #number of files

#make file
for i in $(seq 0 $NUM); do
	echo "Rscript /scripts/4.Count_haplotype_table.R -f hap_data/chr${CHR}/chr${CHR}_haps_${i}.hap -w 5 -c ${CHR} -v hap_data/chr${CHR}/chr${CHR}_haps_${i}.snps -o hap_data/chr${CHR}/chr${CHR}_haps_${i}" >> do_hapcounts_chr${CHR}
done
