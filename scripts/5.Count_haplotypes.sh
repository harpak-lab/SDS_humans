## Count UKB haplotypes (from .haps files)
## updated 5/7/24
## J.M. Cole

##############################################

# Set directories (local cluster)
HAPS=/scratch/ukb/data/		 		# haplotype files 
SCRIPTS=/scratch/ukb/scripts/ 	 	# scripts

#######################
cd /scratch/ukb/data

#make directories
for i in $(seq 1 22); do
	mkdir chr{$i}
done

#move data to folders
for i in $(seq 1 22); do
	mv chr${i}_haps.haps chr${i}
	mv chr${i}_haps_lrs.sample chr${i}
	mv chr${i}.snps.txt chr${i}
done

#Split haplotype files into chunks (150-bp)
for i in $(seq 1 22); do
 echo "split -l 150 -a 3 -d $HAPS/chr${i}/chr${i}_haps.haps $HAPS/chr${i}/chr${i}_haps_" >> split_chroms
done

parallel --jobs 70 < split_chroms

## Format haplotype file names
#Get rid of leading zeros
for i in {1..22}; do
    cd chr$i
    for FILE in *; do 
        mv "$FILE" "$(echo $FILE | sed -e 's/_0*/_/g')"
    done
done

#restore the zero on the first file
for i in $(seq 1 22); do
	mv chr${i}/chr${i}_haps_ chr${i}/chr${i}_haps_0
done


######################################################################

#### Transpose haplotype files
bash $SCRIPTS/Transpose_haplotypes.sh 1 289
bash $SCRIPTS/Transpose_haplotypes.sh 2 294
bash $SCRIPTS/Transpose_haplotypes.sh 3 247
bash $SCRIPTS/Transpose_haplotypes.sh 4 232
bash $SCRIPTS/Transpose_haplotypes.sh 5 217
bash $SCRIPTS/Transpose_haplotypes.sh 6 253
bash $SCRIPTS/Transpose_haplotypes.sh 7 201
bash $SCRIPTS/Transpose_haplotypes.sh 8 189
bash $SCRIPTS/Transpose_haplotypes.sh 9 161
bash $SCRIPTS/Transpose_haplotypes.sh 10 183
bash $SCRIPTS/Transpose_haplotypes.sh 11 181
bash $SCRIPTS/Transpose_haplotypes.sh 12 173
bash $SCRIPTS/Transpose_haplotypes.sh 13 126
bash $SCRIPTS/Transpose_haplotypes.sh 14 118
bash $SCRIPTS/Transpose_haplotypes.sh 15 115
bash $SCRIPTS/Transpose_haplotypes.sh 16 128
bash $SCRIPTS/Transpose_haplotypes.sh 17 119
bash $SCRIPTS/Transpose_haplotypes.sh 18 110
bash $SCRIPTS/Transpose_haplotypes.sh 19 95
bash $SCRIPTS/Transpose_haplotypes.sh 20 96
bash $SCRIPTS/Transpose_haplotypes.sh 21 55
bash $SCRIPTS/Transpose_haplotypes.sh 22 59

cat do_transpose_chr{1..22} > do_transpose

parallel --jobs 70 < do_transpose 


############# Count haplotypes

## haplotype counts (observed)
bash $SCRIPTS/jobs_hapcounts.sh 1 289
bash $SCRIPTS/jobs_hapcounts.sh 2 294
bash $SCRIPTS/jobs_hapcounts.sh 3 247
bash $SCRIPTS/jobs_hapcounts.sh 4 232
bash $SCRIPTS/jobs_hapcounts.sh 5 217
bash $SCRIPTS/jobs_hapcounts.sh 6 253
bash $SCRIPTS/jobs_hapcounts.sh 7 201
bash $SCRIPTS/jobs_hapcounts.sh 8 189
bash $SCRIPTS/jobs_hapcounts.sh 9 161
bash $SCRIPTS/jobs_hapcounts.sh 10 183
bash $SCRIPTS/jobs_hapcounts.sh 11 181
bash $SCRIPTS/jobs_hapcounts.sh 12 173
bash $SCRIPTS/jobs_hapcounts.sh 13 126
bash $SCRIPTS/jobs_hapcounts.sh 14 118
bash $SCRIPTS/jobs_hapcounts.sh 15 115
bash $SCRIPTS/jobs_hapcounts.sh 16 128
bash $SCRIPTS/jobs_hapcounts.sh 17 119
bash $SCRIPTS/jobs_hapcounts.sh 18 110
bash $SCRIPTS/jobs_hapcounts.sh 19 95
bash $SCRIPTS/jobs_hapcounts.sh 20 96
bash $SCRIPTS/jobs_hapcounts.sh 21 55
bash $SCRIPTS/jobs_hapcounts.sh 22 59

cat do_hapcounts_chr{1..22} > do_hapcounts

parallel --jobs 60 < do_hapcounts


## haplotype counts (permuted)
bash $SCRIPTS/jobs_permcounts.sh 1 289
bash $SCRIPTS/jobs_permcounts.sh 2 294
bash $SCRIPTS/jobs_permcounts.sh 3 247
bash $SCRIPTS/jobs_permcounts.sh 4 232
bash $SCRIPTS/jobs_permcounts.sh 5 217
bash $SCRIPTS/jobs_permcounts.sh 6 253
bash $SCRIPTS/jobs_permcounts.sh 7 201
bash $SCRIPTS/jobs_permcounts.sh 8 189
bash $SCRIPTS/jobs_permcounts.sh 9 161
bash $SCRIPTS/jobs_permcounts.sh 10 183
bash $SCRIPTS/jobs_permcounts.sh 11 181
bash $SCRIPTS/jobs_permcounts.sh 12 173
bash $SCRIPTS/jobs_permcounts.sh 13 126
bash $SCRIPTS/jobs_permcounts.sh 14 118
bash $SCRIPTS/jobs_permcounts.sh 15 115
bash $SCRIPTS/jobs_permcounts.sh 16 128
bash $SCRIPTS/jobs_permcounts.sh 17 119
bash $SCRIPTS/jobs_permcounts.sh 18 110
bash $SCRIPTS/jobs_permcounts.sh 19 95
bash $SCRIPTS/jobs_permcounts.sh 20 96
bash $SCRIPTS/jobs_permcounts.sh 21 55
bash $SCRIPTS/jobs_permcounts.sh 22 59

cat do_permcounts_chr{1..22} > do_permcounts

parallel --jobs 70 < do_permcounts


### GET TOTAL HAPCOUNTS ACROSS CHROMS
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.counts $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.counts $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.counts
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.counts $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.counts > $HAPS/chr${i}/chr${i}_haps.counts
	rm -r $HAPS/chr${i}/copy/tmp
done

cp $HAPS/chr*/chr*_haps.counts $HAPS/
mkdir $HAPS/tmp
mv $HAPS/chr1_haps.counts $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.counts
mv $HAPS/tmp/chr1_haps.counts $HAPS/
cat $HAPS/chr{1..22}_haps.counts > $HAPS/all_Chroms.counts

### GET TOTAL HAPCOUNTS ACROSS CHROMS (PERMUTED)
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.counts $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.counts $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.counts
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.counts $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.counts > $HAPS/chr${i}/chr${i}_haps.permuted.counts
	rm -r $HAPS/chr${i}/copy/tmp
done

cp $HAPS/chr*/chr*_haps.permuted.counts $HAPS/
mkdir $HAPS/tmp
mv $HAPS/chr1_haps.permuted.counts $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.permuted.counts
mv $HAPS/tmp/chr1_haps.permuted.counts $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.counts > $HAPS/all_Chroms.permuted.counts