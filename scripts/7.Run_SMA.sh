## Run SMA regressions on 27 traits
## updated 5/10/24
## J.M. Cole

##############################################

# Set directories (local cluster)
DATA=/scratch/ukb/data/		 		# path to data
SCRIPTS=/scratch/ukb/scripts/ 	 	# scripts

####################### 

#get datasets for each trait, job files
cat $DATA/traits/trait_list.txt | while read line 
do
	echo "Rscript $SCRIPTS/Combine_GWAS_SAS.R ${line}" >> $DATA/traits/get_trait_data
done

#run
parallel --jobs 27 < get_trait_data

## do SMA regressions, job files
cat $DATA/traits/trait_list.txt | while read line 
do
	echo "Rscript $SCRIPTS/Perform_regression_SMA.R ${line} 1 1" >> $DATA/do_regressions
done

#run
parallel --jobs 27 < do_regressions