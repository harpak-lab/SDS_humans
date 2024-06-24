## Run likelihood analyses and bootstrapping, local cluster
## updated 5/7/24
## J.M. Cole

##############################################

# Set directories (local cluster)
HAPS=/scratch/ukb/data/		 		# haplotype files 
SCRIPTS=/scratch/ukb/scripts/ 	 	# scripts

####################### Create Job files and run

### ML - viability - create job files
bash $SCRIPTS/jobs_doMLvia.sh 1 289
bash $SCRIPTS/jobs_doMLvia.sh 2 294
bash $SCRIPTS/jobs_doMLvia.sh 3 247
bash $SCRIPTS/jobs_doMLvia.sh 4 232
bash $SCRIPTS/jobs_doMLvia.sh 5 217
bash $SCRIPTS/jobs_doMLvia.sh 6 253
bash $SCRIPTS/jobs_doMLvia.sh 7 201
bash $SCRIPTS/jobs_doMLvia.sh 8 189
bash $SCRIPTS/jobs_doMLvia.sh 9 161
bash $SCRIPTS/jobs_doMLvia.sh 10 183
bash $SCRIPTS/jobs_doMLvia.sh 11 181
bash $SCRIPTS/jobs_doMLvia.sh 12 173
bash $SCRIPTS/jobs_doMLvia.sh 13 126
bash $SCRIPTS/jobs_doMLvia.sh 14 118
bash $SCRIPTS/jobs_doMLvia.sh 15 115
bash $SCRIPTS/jobs_doMLvia.sh 16 128
bash $SCRIPTS/jobs_doMLvia.sh 17 119
bash $SCRIPTS/jobs_doMLvia.sh 18 110
bash $SCRIPTS/jobs_doMLvia.sh 19 95
bash $SCRIPTS/jobs_doMLvia.sh 20 96
bash $SCRIPTS/jobs_doMLvia.sh 21 55
bash $SCRIPTS/jobs_doMLvia.sh 22 59

cat do_ml_chr{1..22} > do_ml_all

#run jobs
parallel --jobs 100 < do_ml_all

### ML- viability (permuted) - create job files
bash $SCRIPTS/jobs_doMLvia_perm.sh 1 289
bash $SCRIPTS/jobs_doMLvia_perm.sh 2 294
bash $SCRIPTS/jobs_doMLvia_perm.sh 3 247
bash $SCRIPTS/jobs_doMLvia_perm.sh 4 232
bash $SCRIPTS/jobs_doMLvia_perm.sh 5 217
bash $SCRIPTS/jobs_doMLvia_perm.sh 6 253
bash $SCRIPTS/jobs_doMLvia_perm.sh 7 201
bash $SCRIPTS/jobs_doMLvia_perm.sh 8 189
bash $SCRIPTS/jobs_doMLvia_perm.sh 9 161
bash $SCRIPTS/jobs_doMLvia_perm.sh 10 183
bash $SCRIPTS/jobs_doMLvia_perm.sh 11 181
bash $SCRIPTS/jobs_doMLvia_perm.sh 12 173
bash $SCRIPTS/jobs_doMLvia_perm.sh 13 126
bash $SCRIPTS/jobs_doMLvia_perm.sh 14 118
bash $SCRIPTS/jobs_doMLvia_perm.sh 15 115
bash $SCRIPTS/jobs_doMLvia_perm.sh 16 128
bash $SCRIPTS/jobs_doMLvia_perm.sh 17 119
bash $SCRIPTS/jobs_doMLvia_perm.sh 18 110
bash $SCRIPTS/jobs_doMLvia_perm.sh 19 95
bash $SCRIPTS/jobs_doMLvia_perm.sh 20 96
bash $SCRIPTS/jobs_doMLvia_perm.sh 21 55
bash $SCRIPTS/jobs_doMLvia_perm.sh 22 59

cat do_ml_perm_chr{1..22} > do_ml_perm_all

#run jobs
parallel --jobs 100 < do_ml_perm_all


########## Post-processing 
#obs
for i in $(seq 1 22); do
	mkdir $HAPS/chr${i}/copy
	cp $HAPS/chr${i}/chr${i}_haps_*.v.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.v.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.v.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.v.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.v.maxl > $HAPS/chr${i}/copy/chr${i}_haps.v.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.v.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.v.maxl
done

#perms
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.v.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.v.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.v.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.v.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.v.maxl > $HAPS/chr${i}/copy/chr${i}_haps.permuted.v.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.v.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.permuted.v.maxl
done

cp $HAPS/chr*/chr*_haps.v.maxl.combined $HAPS/
mkdir $HAPS/tmp
mv $HAPS/chr1_haps.v.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.v.maxl.combined
mv $HAPS/tmp/chr1_haps.v.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.v.maxl.combined > $HAPS/all_Chroms.v.maxl.combined

cp $HAPS/chr*/chr*_haps.permuted.v.maxl.combined $HAPS/
mv $HAPS/chr1_haps.permuted.v.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.permuted.v.maxl.combined
mv $HAPS/tmp/chr1_haps.permuted.v.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.v.maxl.combined > $HAPS/all_Chroms.permuted.v.maxl.combined

### if needing to remove 
for i in $(seq 1 22); do
	rm $HAPS/chr${i}/*.combined
	rm $HAPS/chr${i}/*.maxl
	rm -r $HAPS/chr${i}/copy
done

# check for missing files
bash $SCRIPTS/check_missing.sh$HAPS/chr1/chr1_haps_ 289 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr22/chr22_haps_ 59 .v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr1/chr1_haps_ 289 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .permuted.v.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr22_haps_ 59 .permuted.v.maxl

### BOOTSTRAPPING - viability - create job files
bash $SCRIPTS/jobs_dobootvia.sh 1 289
bash $SCRIPTS/jobs_dobootvia.sh 2 294
bash $SCRIPTS/jobs_dobootvia.sh 3 247
bash $SCRIPTS/jobs_dobootvia.sh 4 232
bash $SCRIPTS/jobs_dobootvia.sh 5 217
bash $SCRIPTS/jobs_dobootvia.sh 6 253
bash $SCRIPTS/jobs_dobootvia.sh 7 201
bash $SCRIPTS/jobs_dobootvia.sh 8 189
bash $SCRIPTS/jobs_dobootvia.sh 9 161
bash $SCRIPTS/jobs_dobootvia.sh 10 183
bash $SCRIPTS/jobs_dobootvia.sh 11 181
bash $SCRIPTS/jobs_dobootvia.sh 12 173
bash $SCRIPTS/jobs_dobootvia.sh 13 126
bash $SCRIPTS/jobs_dobootvia.sh 14 118
bash $SCRIPTS/jobs_dobootvia.sh 15 115
bash $SCRIPTS/jobs_dobootvia.sh 16 128
bash $SCRIPTS/jobs_dobootvia.sh 17 119
bash $SCRIPTS/jobs_dobootvia.sh 18 110
bash $SCRIPTS/jobs_dobootvia.sh 19 95
bash $SCRIPTS/jobs_dobootvia.sh 20 96
bash $SCRIPTS/jobs_dobootvia.sh 21 55
bash $SCRIPTS/jobs_dobootvia.sh 22 59

cat do_boot_chr{1..22} > do_boot_all

#run jobs
parallel --jobs 100 < do_boot_all

####### CHECK THAT BOOTSTRAP FILES GO TO 100 PER WINDOW, run in relevant folder
# Loop through all chr*_haps_*.v.boot files
for file in chr*_haps_*.v.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

##### GET VIABILITY RESULTS - create job files
bash $SCRIPTS/jobs_getresult_via.sh 1 289
bash $SCRIPTS/jobs_getresult_via.sh 2 294
bash $SCRIPTS/jobs_getresult_via.sh 3 247
bash $SCRIPTS/jobs_getresult_via.sh 4 232
bash $SCRIPTS/jobs_getresult_via.sh 5 217
bash $SCRIPTS/jobs_getresult_via.sh 6 253
bash $SCRIPTS/jobs_getresult_via.sh 7 201
bash $SCRIPTS/jobs_getresult_via.sh 8 189
bash $SCRIPTS/jobs_getresult_via.sh 9 161
bash $SCRIPTS/jobs_getresult_via.sh 10 183
bash $SCRIPTS/jobs_getresult_via.sh 11 181
bash $SCRIPTS/jobs_getresult_via.sh 12 173
bash $SCRIPTS/jobs_getresult_via.sh 13 126
bash $SCRIPTS/jobs_getresult_via.sh 14 118
bash $SCRIPTS/jobs_getresult_via.sh 15 115
bash $SCRIPTS/jobs_getresult_via.sh 16 128
bash $SCRIPTS/jobs_getresult_via.sh 17 119
bash $SCRIPTS/jobs_getresult_via.sh 18 110
bash $SCRIPTS/jobs_getresult_via.sh 19 95
bash $SCRIPTS/jobs_getresult_via.sh 20 96
bash $SCRIPTS/jobs_getresult_via.sh 21 55
bash $SCRIPTS/jobs_getresult_via.sh 22 59

cat do_getresult_via_chr{1..22} > do_getresult_via_all 

#run jobs
parallel --jobs 100 < do_getresult_via_all 

##Post processing of results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.v.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.v.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.v.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.v.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.v.result > $HAPS/chr${i}/copy/chr${i}_haps.v.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.v.raw.result -o $HAPS/chr${i}/chr${i}_haps.v.result
	mv $HAPS/chr${i}/chr${i}_haps.v.result.combined $HAPS/chr${i}/chr${i}_haps.v.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.v.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.v.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.v.result
mv $HAPS/tmp/chr1_haps.v.result $HAPS/
cat $HAPS/chr{1..22}_haps.v.result > $HAPS/allchroms.v.result
rm -r $HAPS/tmp 

#### Bootstrapping viability - permuted distribution - create job files
bash $SCRIPTS/jobs_dobootvia_perm.sh 1 289
bash $SCRIPTS/jobs_dobootvia_perm.sh 2 294
bash $SCRIPTS/jobs_dobootvia_perm.sh 3 247
bash $SCRIPTS/jobs_dobootvia_perm.sh 4 232
bash $SCRIPTS/jobs_dobootvia_perm.sh 5 217
bash $SCRIPTS/jobs_dobootvia_perm.sh 6 253
bash $SCRIPTS/jobs_dobootvia_perm.sh 7 201
bash $SCRIPTS/jobs_dobootvia_perm.sh 8 189
bash $SCRIPTS/jobs_dobootvia_perm.sh 9 161
bash $SCRIPTS/jobs_dobootvia_perm.sh 10 183
bash $SCRIPTS/jobs_dobootvia_perm.sh 11 181
bash $SCRIPTS/jobs_dobootvia_perm.sh 12 173
bash $SCRIPTS/jobs_dobootvia_perm.sh 13 126
bash $SCRIPTS/jobs_dobootvia_perm.sh 14 118
bash $SCRIPTS/jobs_dobootvia_perm.sh 15 115
bash $SCRIPTS/jobs_dobootvia_perm.sh 16 128
bash $SCRIPTS/jobs_dobootvia_perm.sh 17 119
bash $SCRIPTS/jobs_dobootvia_perm.sh 18 110
bash $SCRIPTS/jobs_dobootvia_perm.sh 19 95
bash $SCRIPTS/jobs_dobootvia_perm.sh 20 96
bash $SCRIPTS/jobs_dobootvia_perm.sh 21 55
bash $SCRIPTS/jobs_dobootvia_perm.sh 22 59

cat do_getresult_via_chr{1..22} > do_getresult_via_all 

#run jobs
parallel --jobs 100 < do_getresult_via_all 

####### CHECK THAT PERMUTED BOOTSTRAP FILES GO TO 100 PER WINDOW, run in relevant folder
# Loop through all chr*_haps_*.v.boot files
for file in chr*_haps_*.permuted.v.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

## GET VIABILITY PERMUTED RESULTS
bash $SCRIPTS/jobs_getresult_perm_via.sh 1 289
bash $SCRIPTS/jobs_getresult_perm_via.sh 2 294
bash $SCRIPTS/jobs_getresult_perm_via.sh 3 247
bash $SCRIPTS/jobs_getresult_perm_via.sh 4 232
bash $SCRIPTS/jobs_getresult_perm_via.sh 5 217
bash $SCRIPTS/jobs_getresult_perm_via.sh 6 253
bash $SCRIPTS/jobs_getresult_perm_via.sh 7 201
bash $SCRIPTS/jobs_getresult_perm_via.sh 8 189
bash $SCRIPTS/jobs_getresult_perm_via.sh 9 161
bash $SCRIPTS/jobs_getresult_perm_via.sh 10 183
bash $SCRIPTS/jobs_getresult_perm_via.sh 11 181
bash $SCRIPTS/jobs_getresult_perm_via.sh 12 173
bash $SCRIPTS/jobs_getresult_perm_via.sh 13 126
bash $SCRIPTS/jobs_getresult_perm_via.sh 14 118
bash $SCRIPTS/jobs_getresult_perm_via.sh 15 115
bash $SCRIPTS/jobs_getresult_perm_via.sh 16 128
bash $SCRIPTS/jobs_getresult_perm_via.sh 17 119
bash $SCRIPTS/jobs_getresult_perm_via.sh 18 110
bash $SCRIPTS/jobs_getresult_perm_via.sh 19 95
bash $SCRIPTS/jobs_getresult_perm_via.sh 20 96
bash $SCRIPTS/jobs_getresult_perm_via.sh 21 55
bash $SCRIPTS/jobs_getresult_perm_via.sh 22 59

cat do_getresult_perm_via_chr{1..22} > do_getresult_perm_via_all

#run jobs
parallel --jobs 100 < do_getresult_perm_via_all

##Post processing of permuted bootstrapped results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.v.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.v.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.v.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.v.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.v.result > $HAPS/chr${i}/copy/chr${i}_haps.permuted.v.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.v.raw.result -o $HAPS/chr${i}/chr${i}_haps.permuted.v.result
	mv $HAPS/chr${i}/chr${i}_haps.permuted.v.result.combined $HAPS/chr${i}/chr${i}_haps.permuted.v.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.permuted.v.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.permuted.v.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.permuted.v.result
mv $HAPS/tmp/chr1_haps.permuted.v.result $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.v.result > $HAPS/allchroms.permuted.v.result
rm -r $HAPS/tmp 

##############################################################################################
## Run total selection

#TOTAL SELECTION ML - create job files
bash $SCRIPTS/jobs_doMLtot.sh 1 289
bash $SCRIPTS/jobs_doMLtot.sh 2 294
bash $SCRIPTS/jobs_doMLtot.sh 3 247
bash $SCRIPTS/jobs_doMLtot.sh 4 232
bash $SCRIPTS/jobs_doMLtot.sh 5 217
bash $SCRIPTS/jobs_doMLtot.sh 6 253
bash $SCRIPTS/jobs_doMLtot.sh 7 201
bash $SCRIPTS/jobs_doMLtot.sh 8 189
bash $SCRIPTS/jobs_doMLtot.sh 9 161
bash $SCRIPTS/jobs_doMLtot.sh 10 183
bash $SCRIPTS/jobs_doMLtot.sh 11 181
bash $SCRIPTS/jobs_doMLtot.sh 12 173
bash $SCRIPTS/jobs_doMLtot.sh 13 126
bash $SCRIPTS/jobs_doMLtot.sh 14 118
bash $SCRIPTS/jobs_doMLtot.sh 15 115
bash $SCRIPTS/jobs_doMLtot.sh 16 128
bash $SCRIPTS/jobs_doMLtot.sh 17 119
bash $SCRIPTS/jobs_doMLtot.sh 18 110
bash $SCRIPTS/jobs_doMLtot.sh 19 95
bash $SCRIPTS/jobs_doMLtot.sh 20 96
bash $SCRIPTS/jobs_doMLtot.sh 21 55
bash $SCRIPTS/jobs_doMLtot.sh 22 59

cat do_ml_t_chr{1..22} > do_ml_t_all

#run jobs
parallel --jobs 100 < do_ml_t_all

#TOTAL SELECTION PERMUTE
bash $SCRIPTS/jobs_doMLtot_perm.sh 1 289
bash $SCRIPTS/jobs_doMLtot_perm.sh 2 294
bash $SCRIPTS/jobs_doMLtot_perm.sh 3 247
bash $SCRIPTS/jobs_doMLtot_perm.sh 4 232
bash $SCRIPTS/jobs_doMLtot_perm.sh 5 217
bash $SCRIPTS/jobs_doMLtot_perm.sh 6 253
bash $SCRIPTS/jobs_doMLtot_perm.sh 7 201
bash $SCRIPTS/jobs_doMLtot_perm.sh 8 189
bash $SCRIPTS/jobs_doMLtot_perm.sh 9 161
bash $SCRIPTS/jobs_doMLtot_perm.sh 10 183
bash $SCRIPTS/jobs_doMLtot_perm.sh 11 181
bash $SCRIPTS/jobs_doMLtot_perm.sh 12 173
bash $SCRIPTS/jobs_doMLtot_perm.sh 13 126
bash $SCRIPTS/jobs_doMLtot_perm.sh 14 118
bash $SCRIPTS/jobs_doMLtot_perm.sh 15 115
bash $SCRIPTS/jobs_doMLtot_perm.sh 16 128
bash $SCRIPTS/jobs_doMLtot_perm.sh 17 119
bash $SCRIPTS/jobs_doMLtot_perm.sh 18 110
bash $SCRIPTS/jobs_doMLtot_perm.sh 19 95
bash $SCRIPTS/jobs_doMLtot_perm.sh 20 96
bash $SCRIPTS/jobs_doMLtot_perm.sh 21 55
bash $SCRIPTS/jobs_doMLtot_perm.sh 22 59

cat do_ml_t_perm_chr{1..22} > do_ml_t_perm_all

#run jobs
parallel --jobs 100 < do_ml_t_perm_all

########## Post-processing after running total selection
for i in $(seq 1 22); do
	mkdir $HAPS/chr${i}/copy
	cp $HAPS/chr${i}/chr${i}_haps_*.t.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.t.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.t.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.t.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.t.maxl > $HAPS/chr${i}/copy/chr${i}_haps.t.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.t.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.t.maxl
done

#tot perms
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.t.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.t.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.t.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.t.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.t.maxl > $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.permuted.t.maxl
done

cp $HAPS/chr*/chr*_haps.t.maxl.combined $HAPS/
mkdir $HAPS/tmp
mv $HAPS/chr1_haps.t.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.t.maxl.combined
mv $HAPS/tmp/chr1_haps.t.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.t.maxl.combined > $HAPS/all_Chroms.t.maxl.combined

cp $HAPS/chr*/chr*_haps.permuted.t.maxl.combined $HAPS/
mv $HAPS/chr1_haps.permuted.t.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.permuted.t.maxl.combined
mv $HAPS/tmp/chr1_haps.permuted.t.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.t.maxl.combined > $HAPS/all_Chroms.permuted.t.maxl.combined

# check for missing files
bash $SCRIPTS/check_missing.sh $HAPS/chr1/chr1_haps_ 289 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr22/chr22_haps_ 59 .t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr1/chr1_haps_ 289 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .permuted.t.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr22/chr22_haps_ 59 .permuted.t.maxl

### BOOTSTRAPPING - total - create job files
bash $SCRIPTS/jobs_doboottot.sh 1 289
bash $SCRIPTS/jobs_doboottot.sh 2 294
bash $SCRIPTS/jobs_doboottot.sh 3 247
bash $SCRIPTS/jobs_doboottot.sh 4 232
bash $SCRIPTS/jobs_doboottot.sh 5 217
bash $SCRIPTS/jobs_doboottot.sh 6 253
bash $SCRIPTS/jobs_doboottot.sh 7 201
bash $SCRIPTS/jobs_doboottot.sh 8 189
bash $SCRIPTS/jobs_doboottot.sh 9 161
bash $SCRIPTS/jobs_doboottot.sh 10 183
bash $SCRIPTS/jobs_doboottot.sh 11 181
bash $SCRIPTS/jobs_doboottot.sh 12 173
bash $SCRIPTS/jobs_doboottot.sh 13 126
bash $SCRIPTS/jobs_doboottot.sh 14 118
bash $SCRIPTS/jobs_doboottot.sh 15 115
bash $SCRIPTS/jobs_doboottot.sh 16 128
bash $SCRIPTS/jobs_doboottot.sh 17 119
bash $SCRIPTS/jobs_doboottot.sh 18 110
bash $SCRIPTS/jobs_doboottot.sh 19 95
bash $SCRIPTS/jobs_doboottot.sh 20 96
bash $SCRIPTS/jobs_doboottot.sh 21 55
bash $SCRIPTS/jobs_doboottot.sh 22 59

cat do_boot_tot_chr{1..22} > do_boot_tot_all

#run jobs
parallel --jobs 100 < do_boot_tot_all

####### CHECK THAT BOOTSTRAP FILES GO TO 100 PER WINDOW, run in relevant folder
# Loop through all chr*_haps_*.t.boot files
for file in chr*_haps_*.t.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

##### GET TOTAL RESULTS - create job files
bash $SCRIPTS/jobs_getresult_tot.sh 1 289
bash $SCRIPTS/jobs_getresult_tot.sh 2 294
bash $SCRIPTS/jobs_getresult_tot.sh 3 247
bash $SCRIPTS/jobs_getresult_tot.sh 4 232
bash $SCRIPTS/jobs_getresult_tot.sh 5 217
bash $SCRIPTS/jobs_getresult_tot.sh 6 253
bash $SCRIPTS/jobs_getresult_tot.sh 7 201
bash $SCRIPTS/jobs_getresult_tot.sh 8 189
bash $SCRIPTS/jobs_getresult_tot.sh 9 161
bash $SCRIPTS/jobs_getresult_tot.sh 10 183
bash $SCRIPTS/jobs_getresult_tot.sh 11 181
bash $SCRIPTS/jobs_getresult_tot.sh 12 173
bash $SCRIPTS/jobs_getresult_tot.sh 13 126
bash $SCRIPTS/jobs_getresult_tot.sh 14 118
bash $SCRIPTS/jobs_getresult_tot.sh 15 115
bash $SCRIPTS/jobs_getresult_tot.sh 16 128
bash $SCRIPTS/jobs_getresult_tot.sh 17 119
bash $SCRIPTS/jobs_getresult_tot.sh 18 110
bash $SCRIPTS/jobs_getresult_tot.sh 19 95
bash $SCRIPTS/jobs_getresult_tot.sh 20 96
bash $SCRIPTS/jobs_getresult_tot.sh 21 55
bash $SCRIPTS/jobs_getresult_tot.sh 22 59

cat do_getresult_tot_chr{1..22} > do_getresult_tot_all

#run jobs
parallel --jobs 100 < do_getresult_tot_all

##Post processing of results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.t.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.t.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.t.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.t.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.t.result > $HAPS/chr${i}/copy/chr${i}_haps.t.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.t.raw.result -o $HAPS/chr${i}/chr${i}_haps.t.result
	mv $HAPS/chr${i}/chr${i}_haps.t.result.combined $HAPS/chr${i}/chr${i}_haps.t.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.t.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.t.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.t.result
mv $HAPS/tmp/chr1_haps.t.result $HAPS/
cat $HAPS/chr{1..22}_haps.t.result > $HAPS/allchroms.t.result
rm -r $HAPS/tmp 

#### Bootstrapping total - permuted distribution - create job files
bash $SCRIPTS/jobs_doboottot_perm.sh 1 289
bash $SCRIPTS/jobs_doboottot_perm.sh 2 294
bash $SCRIPTS/jobs_doboottot_perm.sh 3 247
bash $SCRIPTS/jobs_doboottot_perm.sh 4 232
bash $SCRIPTS/jobs_doboottot_perm.sh 5 217
bash $SCRIPTS/jobs_doboottot_perm.sh 6 253
bash $SCRIPTS/jobs_doboottot_perm.sh 7 201
bash $SCRIPTS/jobs_doboottot_perm.sh 8 189
bash $SCRIPTS/jobs_doboottot_perm.sh 9 161
bash $SCRIPTS/jobs_doboottot_perm.sh 10 183
bash $SCRIPTS/jobs_doboottot_perm.sh 11 181
bash $SCRIPTS/jobs_doboottot_perm.sh 12 173
bash $SCRIPTS/jobs_doboottot_perm.sh 13 126
bash $SCRIPTS/jobs_doboottot_perm.sh 14 118
bash $SCRIPTS/jobs_doboottot_perm.sh 15 115
bash $SCRIPTS/jobs_doboottot_perm.sh 16 128
bash $SCRIPTS/jobs_doboottot_perm.sh 17 119
bash $SCRIPTS/jobs_doboottot_perm.sh 18 110
bash $SCRIPTS/jobs_doboottot_perm.sh 19 95
bash $SCRIPTS/jobs_doboottot_perm.sh 20 96
bash $SCRIPTS/jobs_doboottot_perm.sh 21 55
bash $SCRIPTS/jobs_doboottot_perm.sh 22 59

cat do_permboot_tot_chr{1..22} > do_permboot_tot_all

#run jobs
parallel --jobs 100 < do_permboot_tot_all

# Loop through all chr*_haps_*.permuted.t.boot files
for file in chr*_haps_*.permuted.t.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

## GET TOTAL PERMUTED RESULTS
bash $SCRIPTs/jobs_getresult_perm_tot.sh 1 289
bash $SCRIPTs/jobs_getresult_perm_tot.sh 2 294
bash $SCRIPTs/jobs_getresult_perm_tot.sh 3 247
bash $SCRIPTs/jobs_getresult_perm_tot.sh 4 232
bash $SCRIPTs/jobs_getresult_perm_tot.sh 5 217
bash $SCRIPTs/jobs_getresult_perm_tot.sh 6 253
bash $SCRIPTs/jobs_getresult_perm_tot.sh 7 201
bash $SCRIPTs/jobs_getresult_perm_tot.sh 8 189
bash $SCRIPTs/jobs_getresult_perm_tot.sh 9 161
bash $SCRIPTs/jobs_getresult_perm_tot.sh 10 183
bash $SCRIPTs/jobs_getresult_perm_tot.sh 11 181
bash $SCRIPTs/jobs_getresult_perm_tot.sh 12 173
bash $SCRIPTs/jobs_getresult_perm_tot.sh 13 126
bash $SCRIPTs/jobs_getresult_perm_tot.sh 14 118
bash $SCRIPTs/jobs_getresult_perm_tot.sh 15 115
bash $SCRIPTs/jobs_getresult_perm_tot.sh 16 128
bash $SCRIPTs/jobs_getresult_perm_tot.sh 17 119
bash $SCRIPTs/jobs_getresult_perm_tot.sh 18 110
bash $SCRIPTs/jobs_getresult_perm_tot.sh 19 95
bash $SCRIPTs/jobs_getresult_perm_tot.sh 20 96
bash $SCRIPTs/jobs_getresult_perm_tot.sh 21 55
bash $SCRIPTs/jobs_getresult_perm_tot.sh 22 59

cat do_getresult_perm_tot_chr{1..22} > do_getresult_perm_tot_all

#run jobs
parallel --jobs 100 < do_getresult_perm_tot_all

##Post processing of permuted bootstrapped results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.t.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.t.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.t.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.t.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.t.result > $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.result -o $HAPS/chr${i}/chr${i}_haps.permuted.t.result
	mv $HAPS/chr${i}/chr${i}_haps.permuted.t.result.combined $HAPS/chr${i}/chr${i}_haps.permuted.t.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.permuted.t.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.permuted.t.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.permuted.t.result
mv $HAPS/tmp/chr1_haps.permuted.t.result $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.t.result > $HAPS/allchroms.permuted.t.result
rm -r $HAPS/tmp 


##############################################################################################
## Run fecundity 

#FECUNDITY SELECTION ML - create job files
bash $SCRIPTS/jobs_doMLfec.sh 1 289
bash $SCRIPTS/jobs_doMLfec.sh 2 294
bash $SCRIPTS/jobs_doMLfec.sh 3 247
bash $SCRIPTS/jobs_doMLfec.sh 4 232
bash $SCRIPTS/jobs_doMLfec.sh 5 217
bash $SCRIPTS/jobs_doMLfec.sh 6 253
bash $SCRIPTS/jobs_doMLfec.sh 7 201
bash $SCRIPTS/jobs_doMLfec.sh 8 189
bash $SCRIPTS/jobs_doMLfec.sh 9 161
bash $SCRIPTS/jobs_doMLfec.sh 10 183
bash $SCRIPTS/jobs_doMLfec.sh 11 181
bash $SCRIPTS/jobs_doMLfec.sh 12 173
bash $SCRIPTS/jobs_doMLfec.sh 13 126
bash $SCRIPTS/jobs_doMLfec.sh 14 118
bash $SCRIPTS/jobs_doMLfec.sh 15 115
bash $SCRIPTS/jobs_doMLfec.sh 16 128
bash $SCRIPTS/jobs_doMLfec.sh 17 119
bash $SCRIPTS/jobs_doMLfec.sh 18 110
bash $SCRIPTS/jobs_doMLfec.sh 19 95
bash $SCRIPTS/jobs_doMLfec.sh 20 96
bash $SCRIPTS/jobs_doMLfec.sh 21 55
bash $SCRIPTS/jobs_doMLfec.sh 22 59

cat do_ml_f_chr{1..22} > do_ml_f_all

#run jobs
parallel --jobs 100 < do_ml_f_all


#FECUNDITY SELECTION PERMUTE - create job files
bash $SCRIPTS/jobs_doMLfec_perm.sh 1 289
bash $SCRIPTS/jobs_doMLfec_perm.sh 2 294
bash $SCRIPTS/jobs_doMLfec_perm.sh 3 247
bash $SCRIPTS/jobs_doMLfec_perm.sh 4 232
bash $SCRIPTS/jobs_doMLfec_perm.sh 5 217
bash $SCRIPTS/jobs_doMLfec_perm.sh 6 253
bash $SCRIPTS/jobs_doMLfec_perm.sh 7 201
bash $SCRIPTS/jobs_doMLfec_perm.sh 8 189
bash $SCRIPTS/jobs_doMLfec_perm.sh 9 161
bash $SCRIPTS/jobs_doMLfec_perm.sh 10 183
bash $SCRIPTS/jobs_doMLfec_perm.sh 11 181
bash $SCRIPTS/jobs_doMLfec_perm.sh 12 173
bash $SCRIPTS/jobs_doMLfec_perm.sh 13 126
bash $SCRIPTS/jobs_doMLfec_perm.sh 14 118
bash $SCRIPTS/jobs_doMLfec_perm.sh 15 115
bash $SCRIPTS/jobs_doMLfec_perm.sh 16 128
bash $SCRIPTS/jobs_doMLfec_perm.sh 17 119
bash $SCRIPTS/jobs_doMLfec_perm.sh 18 110
bash $SCRIPTS/jobs_doMLfec_perm.sh 19 95
bash $SCRIPTS/jobs_doMLfec_perm.sh 20 96
bash $SCRIPTS/jobs_doMLfec_perm.sh 21 55
bash $SCRIPTS/jobs_doMLfec_perm.sh 22 59

cat do_ml_f_perm_chr{1..22} > do_ml_f_perm_all

#run jobs
parallel --jobs 100 < do_ml_f_perm_all

########## Post-processing after running fecundity selection
for i in $(seq 1 22); do
	mkdir $HAPS/chr${i}/copy
	cp $HAPS/chr${i}/chr${i}_haps_*.f.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.f.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.f.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.f.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.f.maxl > $HAPS/chr${i}/copy/chr${i}_haps.t.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.t.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.f.maxl
done

#fec perms
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.f.maxl $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.f.maxl $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.f.maxl
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.f.maxl $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.f.maxl > $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.maxl
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.t.raw.maxl -o $HAPS/chr${i}/chr${i}_haps.permuted.f.maxl
done

cp $HAPS/chr*/chr*_haps.f.maxl.combined $HAPS/
mkdir $HAPS/tmp
mv $HAPS/chr1_haps.f.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.f.maxl.combined
mv $HAPS/tmp/chr1_haps.f.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.f.maxl.combined > $HAPS/all_Chroms.f.maxl.combined

cp $HAPS/chr*/chr*_haps.permuted.f.maxl.combined $HAPS/
mv $HAPS/chr1_haps.permuted.f.maxl.combined $HAPS/tmp
sed -i '1d' $HAPS/chr*_haps.permuted.f.maxl.combined
mv $HAPS/tmp/chr1_haps.permuted.f.maxl.combined $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.f.maxl.combined > $HAPS/all_Chroms.permuted.f.maxl.combined

# check for missing files
bash $SCRIPTS/check_missing.sh $HAPS/chr1/chr1_haps_ 289 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr22/chr22_haps_ 59 .f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr1/chr1_haps_ 289 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr2/chr2_haps_ 294 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr3/chr3_haps_ 247 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr4/chr4_haps_ 232 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr5/chr5_haps_ 217 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr6/chr6_haps_ 253 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr7/chr7_haps_ 201 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr8/chr8_haps_ 189 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr9/chr9_haps_ 161 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr10/chr10_haps_ 183 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr11/chr11_haps_ 181 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr12/chr12_haps_ 173 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr13/chr13_haps_ 126 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr14/chr14_haps_ 118 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr15/chr15_haps_ 115 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr16/chr16_haps_ 128 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr17/chr17_haps_ 119 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr18/chr18_haps_ 110 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr19/chr19_haps_ 95 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr20/chr20_haps_ 96 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr21/chr21_haps_ 55 .permuted.f.maxl
bash $SCRIPTS/check_missing.sh $HAPS/chr22/chr22_haps_ 59 .permuted.f.maxl

### BOOTSTRAPPING - fecundity - create job files
bash $SCRIPTS/jobs_dobootfec.sh 1 289
bash $SCRIPTS/jobs_dobootfec.sh 2 294
bash $SCRIPTS/jobs_dobootfec.sh 3 247
bash $SCRIPTS/jobs_dobootfec.sh 4 232
bash $SCRIPTS/jobs_dobootfec.sh 5 217
bash $SCRIPTS/jobs_dobootfec.sh 6 253
bash $SCRIPTS/jobs_dobootfec.sh 7 201
bash $SCRIPTS/jobs_dobootfec.sh 8 189
bash $SCRIPTS/jobs_dobootfec.sh 9 161
bash $SCRIPTS/jobs_dobootfec.sh 10 183
bash $SCRIPTS/jobs_dobootfec.sh 11 181
bash $SCRIPTS/jobs_dobootfec.sh 12 173
bash $SCRIPTS/jobs_dobootfec.sh 13 126
bash $SCRIPTS/jobs_dobootfec.sh 14 118
bash $SCRIPTS/jobs_dobootfec.sh 15 115
bash $SCRIPTS/jobs_dobootfec.sh 16 128
bash $SCRIPTS/jobs_dobootfec.sh 17 119
bash $SCRIPTS/jobs_dobootfec.sh 18 110
bash $SCRIPTS/jobs_dobootfec.sh 19 95
bash $SCRIPTS/jobs_dobootfec.sh 20 96
bash $SCRIPTS/jobs_dobootfec.sh 21 55
bash $SCRIPTS/jobs_dobootfec.sh 22 59

cat do_boot_fec_chr{1..22} > do_boot_fec_all

#run jobs
parallel --jobs 100 < do_boot_fec_all

####### CHECK THAT BOOTSTRAP FILES GO TO 100 PER WINDOW, run in relevant folder
# Loop through all chr*_haps_*.f.boot files
for file in chr*_haps_*.f.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

##### GET FECUNDITY RESULTS - create job files
bash $SCRIPTS/jobs_getresult_fec.sh 1 289
bash $SCRIPTS/jobs_getresult_fec.sh 2 294
bash $SCRIPTS/jobs_getresult_fec.sh 3 247
bash $SCRIPTS/jobs_getresult_fec.sh 4 232
bash $SCRIPTS/jobs_getresult_fec.sh 5 217
bash $SCRIPTS/jobs_getresult_fec.sh 6 253
bash $SCRIPTS/jobs_getresult_fec.sh 7 201
bash $SCRIPTS/jobs_getresult_fec.sh 8 189
bash $SCRIPTS/jobs_getresult_fec.sh 9 161
bash $SCRIPTS/jobs_getresult_fec.sh 10 183
bash $SCRIPTS/jobs_getresult_fec.sh 11 181
bash $SCRIPTS/jobs_getresult_fec.sh 12 173
bash $SCRIPTS/jobs_getresult_fec.sh 13 126
bash $SCRIPTS/jobs_getresult_fec.sh 14 118
bash $SCRIPTS/jobs_getresult_fec.sh 15 115
bash $SCRIPTS/jobs_getresult_fec.sh 16 128
bash $SCRIPTS/jobs_getresult_fec.sh 17 119
bash $SCRIPTS/jobs_getresult_fec.sh 18 110
bash $SCRIPTS/jobs_getresult_fec.sh 19 95
bash $SCRIPTS/jobs_getresult_fec.sh 20 96
bash $SCRIPTS/jobs_getresult_fec.sh 21 55
bash $SCRIPTS/jobs_getresult_fec.sh 22 59

cat do_getresult_fec_chr{1..22} > do_getresult_fec_all

#run jobs
parallel --jobs 100 < do_getresult_fec_all

##Post processing of results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.f.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.f.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.f.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.f.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.f.result > $HAPS/chr${i}/copy/chr${i}_haps.f.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.f.raw.result -o $HAPS/chr${i}/chr${i}_haps.f.result
	mv $HAPS/chr${i}/chr${i}_haps.f.result.combined $HAPS/chr${i}/chr${i}_haps.f.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.f.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.f.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.f.result
mv $HAPS/tmp/chr1_haps.f.result $HAPS/
cat $HAPS/chr{1..22}_haps.f.result > $HAPS/allchroms.f.result
rm -r $HAPS/tmp 

#### Bootstrapping fecundity - permuted distribution - create job files
bash $SCRIPTS/jobs_dobootfec_perm.sh 1 289
bash $SCRIPTS/jobs_dobootfec_perm.sh 2 294
bash $SCRIPTS/jobs_dobootfec_perm.sh 3 247
bash $SCRIPTS/jobs_dobootfec_perm.sh 4 232
bash $SCRIPTS/jobs_dobootfec_perm.sh 5 217
bash $SCRIPTS/jobs_dobootfec_perm.sh 6 253
bash $SCRIPTS/jobs_dobootfec_perm.sh 7 201
bash $SCRIPTS/jobs_dobootfec_perm.sh 8 189
bash $SCRIPTS/jobs_dobootfec_perm.sh 9 161
bash $SCRIPTS/jobs_dobootfec_perm.sh 10 183
bash $SCRIPTS/jobs_dobootfec_perm.sh 11 181
bash $SCRIPTS/jobs_dobootfec_perm.sh 12 173
bash $SCRIPTS/jobs_dobootfec_perm.sh 13 126
bash $SCRIPTS/jobs_dobootfec_perm.sh 14 118
bash $SCRIPTS/jobs_dobootfec_perm.sh 15 115
bash $SCRIPTS/jobs_dobootfec_perm.sh 16 128
bash $SCRIPTS/jobs_dobootfec_perm.sh 17 119
bash $SCRIPTS/jobs_dobootfec_perm.sh 18 110
bash $SCRIPTS/jobs_dobootfec_perm.sh 19 95
bash $SCRIPTS/jobs_dobootfec_perm.sh 20 96
bash $SCRIPTS/jobs_dobootfec_perm.sh 21 55
bash $SCRIPTS/jobs_dobootfec_perm.sh 22 59

cat do_permboot_fec_chr{1..22} > do_permboot_fec_all

#run jobs
parallel --jobs 100 < do_permboot_fec_all

####### CHECK THAT PERMUTED BOOTSTRAP FILES GO TO 100 PER WINDOW, run in relevant folder
# Loop through all chr*_haps_*.f.boot files
for file in chr*_haps_*.permuted.f.boot; do
    echo "Checking $file"
    # Assume the file is OK initially
    file_ok=true

    # Get the list of unique windows
    windows=$(awk 'NR>1 {print $1}' "$file" | sort -nu)

    # For each window, check if the maximum bootrep is 100
    for window in $windows; do
        max_bootrep=$(awk -v win="$window" '$1==win {print $NF}' "$file" | sort -nu | tail -n 1)
        if [ "$max_bootrep" -ne 100 ]; then
            # If any window doesn't go up to 100, mark the file as not OK and break
            file_ok=false
            break
        fi
    done

    # Report the file if it's not OK
    if [ "$file_ok" = false ]; then
        echo "$file does not have all windows going up to 100."
    fi
done

## GET FECUNDITY PERMUTED RESULTS - create job  files
bash $SCRIPTS/jobs_getresult_perm_fec.sh 1 289
bash $SCRIPTS/jobs_getresult_perm_fec.sh 2 294
bash $SCRIPTS/jobs_getresult_perm_fec.sh 3 247
bash $SCRIPTS/jobs_getresult_perm_fec.sh 4 232
bash $SCRIPTS/jobs_getresult_perm_fec.sh 5 217
bash $SCRIPTS/jobs_getresult_perm_fec.sh 6 253
bash $SCRIPTS/jobs_getresult_perm_fec.sh 7 201
bash $SCRIPTS/jobs_getresult_perm_fec.sh 8 189
bash $SCRIPTS/jobs_getresult_perm_fec.sh 9 161
bash $SCRIPTS/jobs_getresult_perm_fec.sh 10 183
bash $SCRIPTS/jobs_getresult_perm_fec.sh 11 181
bash $SCRIPTS/jobs_getresult_perm_fec.sh 12 173
bash $SCRIPTS/jobs_getresult_perm_fec.sh 13 126
bash $SCRIPTS/jobs_getresult_perm_fec.sh 14 118
bash $SCRIPTS/jobs_getresult_perm_fec.sh 15 115
bash $SCRIPTS/jobs_getresult_perm_fec.sh 16 128
bash $SCRIPTS/jobs_getresult_perm_fec.sh 17 119
bash $SCRIPTS/jobs_getresult_perm_fec.sh 18 110
bash $SCRIPTS/jobs_getresult_perm_fec.sh 19 95
bash $SCRIPTS/jobs_getresult_perm_fec.sh 20 96
bash $SCRIPTS/jobs_getresult_perm_fec.sh 21 55
bash $SCRIPTS/jobs_getresult_perm_fec.sh 22 59

cat do_getresult_perm_fec_chr{1..22} > do_getresult_perm_fec_all

#run jobs
parallel --jobs 50 < do_getresult_perm_fec_all

##Post processing of permuted bootstrapped results files
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps_*.permuted.f.result $HAPS/chr${i}/copy
	mkdir $HAPS/chr${i}/copy/tmp
	mv $HAPS/chr${i}/copy/chr${i}_haps_0.permuted.f.result $HAPS/chr${i}/copy/tmp
	sed -i '1d' $HAPS/chr${i}/copy/chr${i}_haps_*.permuted.f.result
	mv $HAPS/chr${i}/copy/tmp/chr${i}_haps_0.permuted.f.result $HAPS/chr${i}/copy
	cat $HAPS/chr${i}/copy/chr${i}_haps_{0..300}.permuted.f.result > $HAPS/chr${i}/copy/chr${i}_haps.permuted.f.raw.result
	rm -r $HAPS/chr${i}/copy/tmp
	Rscript $SCRIPTS/fix_window.R -f $HAPS/chr${i}/copy/chr${i}_haps.permuted.f.raw.result -o $HAPS/chr${i}/chr${i}_haps.permuted.f.result
	mv $HAPS/chr${i}/chr${i}_haps.permuted.f.result.combined $HAPS/chr${i}/chr${i}_haps.permuted.f.result
done

#concat results
for i in $(seq 1 22); do
	cp $HAPS/chr${i}/chr${i}_haps.permuted.f.result $HAPS/
done

mkdir $HAPS/tmp 
mv $HAPS/chr1_haps.permuted.f.result $HAPS/tmp 
sed -i '1d' $HAPS/chr*_haps.permuted.f.result
mv $HAPS/tmp/chr1_haps.permuted.f.result $HAPS/
cat $HAPS/chr{1..22}_haps.permuted.f.result > $HAPS/allchroms.permuted.f.result
rm -r $HAPS/tmp 
