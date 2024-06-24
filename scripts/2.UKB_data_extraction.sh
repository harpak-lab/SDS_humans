## 2. DATA EXTRACTION FOR UK Biobank

## Data ukb45020
##############################################

# Decrypt and unpack using provided handlers
ukbmd5 ukb45020.enc
ukb_unpack ukb45020.enc k61666r45020.key

# Full metadata as csv, r compatible outputs, and data dictionary html files
ukb_conv ukb45020.enc_ukb csv
ukb_conv ukb45020.enc_ukb r
ukb_conv ukb45020.enc_ukb docs

# Retrieve autosomal (chrs 1-22) data for phased haplotypes, genotyped SNPs, and imputed SNPs
for i in $(seq 1 22); do
	ukbgene cal -c${i} -ak61666r45020.key
	ukbgene cal -c${i} -m -ak61666r45020.key
	ukbgene hap -c${i} -ak61666r45020.key
	ukbgene hap -c${i} -m -ak61666r45020.key
	ukbgene imp -c${i} -ak61666r45020.key
	ukbgene imp -c${i} -m -ak61666r45020.key
done

# Metadata and relevant fields for filtering and analyses
## Field 21000 ("Ethnic background")
## Field 21022 ("Age at recruitment")
## Field 22006 ("Genetic ethnic grouping")
## Field 22020 ("Used in genetic principal components")
## Field 22001 ("Genetic sex")
## Field 31 ("Sex")
## Field 22027 ("Outliers for missingness and heterozygosity")
## Field 22019 ("Sex chromosome aneuploidy")
## Field 2734 ("Number of live births")
## Field 2405 ("Number of children fathered")

echo "21000" >> metafields1
echo "21022" >> metafields1
echo "22006" >> metafields1
echo "22020" >> metafields1
echo "22001" >> metafields1
echo "31" >> metafields1
echo "22027" >> metafields1
echo "22019" >> metafields1
echo "2734" >> metafields1
echo "2405" >> metafields1

# Extract metadata (9 fields)
ukbconv ukb45020.enc_ukb csv -imetafields1 -oUKB_metadata_subset_1-24