## 4. PLINK 1 and 2 filtering procedures and processing for UK Biobank data
## updated 5/1/24
## J.M. Cole

##############################################

# Set directories 

HAPS=/scratch/06809/jmcole/ukb/data/		 # haplotype files 
SCRIPTS=/scratch/06809/jmcole/ukb/scripts  	 # processing scripts
META=/scratch/06809/jmcole/ukb/data/metadata # metadata
OUT=/scratch/06809/jmcole/ukb/analysis/ 	 # analysis files and outputs

#######################

# Creates job files using TACC's launcher_creator.py script, submits sbatch for each autosome

# Make directories for 22 autosomes
for i in $(seq 1 22); do
   mkdir chr${i}
done

# Extract autosomal data for samples after sample-level QC (see R scripts for sample-level filtering)
for i in $(seq 1 22); do
	echo "plink2 --bgen $HAPS/ukb_hap_chr${i}_v2.bgen ref-first snpid-chr --allow-extra-chr --sample $HAPS/ukb61666_hap_chr${i}_v2_s487220.sample --keep $META/UKB_filtered_samples_list_1-24.txt --export bgen-1.2 --out $OUT/chr${i}/chr${i}_filtered_step1" > sample_filter_chr${i}
	launcher_creator.py -j sample_filter_chr${i} -n sample_filter_chr${i} -q normal -N 1 -t 48:00:00
	sbatch -A Recombining-sex-chro sample_filter_chr${i}.slurm
done


######### QC Marker filtering steps ##############

## Basic QC
### maf > 0.01
### geno at least 95%
### biallelic snps only
### HWE tests p < 10^-6
### remove private SNPs, sex chromosome homologous snps (from Galichon et al, 2012 and Kasimatis et al, 2021)
### see R script for marker QC for pre-processing of sex-chromosome homologous SNPs

for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step1.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step1.sample --maf 0.01 --geno 0.05 --max-alleles 2 --snps-only --hwe 0.000001 --exclude $OUT/excluded_snps/galichon_snps.txt $OUT/excluded_snps/Kasimatis_removed_snps_BLAT.txt $OUT/excluded_snps/bycroft_correspond_snps.txt --export bgen-1.2 --out $OUT/chr${i}/chr${i}_filtered_step2" > qc_filter_chr${i}
	launcher_creator.py -j qc_filter_chr${i} -n qc_filter_chr${i} -q normal -N 1 -t 48:00:00
	sbatch -A Recombining-sex-chro qc_filter_chr${i}.slurm
done

# Get alleles private in males or females
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step2.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step1.sample  --keep-males --max-maf 0 --write-snplist --out $OUT/chr${i}/snp_list_males_chr${i}" > keep_males_${i}
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step2.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step1.sample --keep-females --max-maf 0 --write-snplist --out $OUT/chr${i}/snp_list_females_chr${i}" > keep_females_${i} 
	launcher_creator.py -j keep_males_${i} -n keep_males_${i} -q normal -N 1 -t 24:00:00
	launcher_creator.py -j keep_females_${i} -n keep_females_${i} -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro keep_males_${i}.slurm
	sbatch -A Recombining-sex-chro keep_females_${i}.slurm  
done

# Remove private alleles in males or females
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step2.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step1.sample --exclude $OUT/chr${i}/snp_list_females_chr${i}.snplist $OUT/chr${i}/snp_list_males_chr${i}.snplist --export bgen-1.2 --out $OUT/chr${i}/chr${i}_filtered_step3" > remove_private_${i}
	launcher_creator.py -j remove_private_${i} -n remove_private_${i} -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro remove_private_${i}.slurm
done

#Note: after QC no private alleles remain in the dataset for either sex. 


## Extract SNP and genotype counts for artefact filtering ###
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step3.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step3.sample  --keep-females --geno-counts --out $OUT/chr${i}/chr${i}.females" > gcount_females_${i}
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step3.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step3.sample  --keep-males --geno-counts --out $OUT/chr${i}/chr${i}.males" > gcount_males_${i}
	launcher_creator.py -j gcount_females_${i} -n gcount_females_${i} -q normal -N 1 -t 24:00:00
	launcher_creator.py -j gcount_females_${i} -n gcount_females_${i} -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro gcount_females_${i}.slurm
	sbatch -A Recombining-sex-chro gcount_males_${i}.slurm
done

#move to output folder
scp $OUT/chr*/*.gcount $OUT

#remove headers on files for chromosomes 2-22
for i in $(seq 2 22); do
 sed -i '1d' chr${i}.males.gcount
 sed -i '1d' chr${i}.females.gcount
done

#combine genotype counts
cat chr{1..22}.males.gcount> allchroms.males.gcount
cat chr{1..22}.females.gcount > allchroms.females.gcount

## Remove SNPs flagged as potential artefacts
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step3.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step3.sample --exclude $OUT/genotype_filters_snps_toremove_5-23.txt --export bgen-1.2 --out $OUT/chr${i}/chr${i}_filtered_step4" > remove_qc2_${i}
	launcher_creator.py -j remove_qc2_${i} -n remove_qc2_${i} -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro remove_qc2_${i}.slurm
done

## Obtain allele frequencies for each chromosome
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step4.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --freq --out $OUT/chr${i}/chr${i}.haps.allele_freqs" > calc_allele_freqs_chr${i}
	launcher_creator.py -j calc_allele_freqs_chr${i} -n calc_allele_freqs_chr${i} -q normal -N 1 -t 48:00:00
	sbatch -A Recombining-sex-chro calc_allele_freqs_chr${i}.slurm
done


## Ensure allele frequencies are minor allele frequencies
for i in $(seq 1 22); do 

awk 'BEGIN {
     IFS = OFS = "\t"
  } 
  NR==1{print;next} {
     for (column = 5; column <= NF; ++column) {
        if ($column > 0.5) {
            $column = 1-$5
        }
     }    
     print 
  }         
' $OUT/chr${i}/chr${i}.haps.allele_freqs.afreq > $OUT/chr${i}/chr${i}.mafs.freq

done

## Output pvar for each chromosome
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step4.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --make-just-pvar --out $OUT/chr${i}/chr${i}" > chr${i}_make_pvar
	launcher_creator.py -j chr${i}_make_pvar -n chr${i}_make_pvar -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro chr${i}_make_pvar.slurm
done

#Download gene annotation file for GCRh37 from Ensembl
wget https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz 

#Get gene coordinates for all chromosomes
awk '($4=="gene") && ($5 ~ /protein_coding/) {OFS="\t"; print $1,$2,$3,$4,$6,$7}' Homo_sapiens.GRCh37.87.gtf > GCRh37.genes.bed

#per chromosome coordinates
for i in $(seq 1 22); do
	cat GCRh37.genes.bed | awk -v c=$i '($1 == $c){OFS="\t"; print $1, $2, $3, $4, $5, $6, $7}' > $OUT/chr$[i}/GCRh37.genes.chr${i}.bed
done

#export all imputed allele frequencies
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step4.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --make-just-pvar --out $OUT/chr${i}/chr${i}" > chr${i}_make_pvar
	launcher_creator.py -j chr${i}_make_pvar -n chr${i}_make_pvar -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro chr${i}_make_pvar.slurm
done

## extract r2 between adjacent SNPs
for i in $(seq 1 22); do 
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step4.bgen ref-first --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --r2 --out $OUT/chr${i}/ukb_ld_chr${i}" > chr${i}_ld
	launcher_creator.py -j chr${i}_ld -n chr${i}_ld -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro chr${i}_ld.slurm
done

#Extract filtered haplotype data as 0/1 for each chromosome
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}/chr${i}_filtered_step4.bgen ref-first --oxford-single-chr $i --allow-extra-chr --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --export haps --out $OUT/chr${i}/chr${i}_haps" > chr${i}_haps
	launcher_creator.py -j chr${i}_haps -n chr${i}_haps -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro chr${i}_haps.slurm
done

#Get sample lists with offspring data
for i in $(seq 1 22); do
	echo "Rscript $SCRIPTS/Merge_sample_data.R -f $META/UKB_sample_LRS_1-24.txt -s $OUT/chr${i}/chr${i}_haps.sample -o  $OUT/chr${i}/chr${i}_haps_lrs.sample" > chr${i}_merge_samples
	launcher_creator.py -j chr${i}_merge_samples -n chr${i}_merge_samples -q normal -N 1 -t 24:00:00
	sbatch -A Recombining-sex-chro chr${i}_merge_samples.slurm
done

#LD-prune SNPs
for i in $(seq 1 22); do
	echo "plink2 --bgen $OUT/chr${i}/chr${i}_filtered_step4.bgen ref-first --allow-extra-chr --oxford-single-chr ${i} --sample $OUT/chr${i}/chr${i}_filtered_step4.sample --indep-pairwise 50 5 0.2 --out $OUT/chr${i}/chr${i}_pruned_02" > chr${i}_pruned
		launcher_creator.py -j chr${i}_pruned -n chr${i}_pruned -q normal -N 1 -t 24:00:00
		sbatch -A Recombining-sex-chro chr${i}_pruned.slurm
done

for i in $(seq 1 22); do
	cp $OUT/chr${i}/chr${i}_pruned_02.prune.in $OUT
done

cat $OUT/chr{1..22}_pruned_02.prune.in > $OUT/all_chrs_pruned_02.txt  


## Transfer files to local cluster for analyses
scp $OUT/chr*/chr*_haps.haps $OUT/chr*/chr*_haps_lrs.sample $OUT/chr*/chr*.snps.txt $OUT/all_chrs_pruned_02.txt jmcole@ccbb:/scratch/ukb/data/
