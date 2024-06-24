## SITE_LEVEL FILTERING OF UK BIOBANK
## updated 5/17/24
## J.M. Cole
##
## Identify sex chromosome homologous SNPs
## Tests of between-sex missingingness, inflated heterozygosity, and deficits of minar allele homozygosity
## See Ruzicka et al, 2022, PloS Biol
###########################################################################
setwd("/home/jared/Documents/Projects/SAS/humans/UKB/downloaded_data/")

### Homology detected via BLAT in Kasimatis et al, 2021 (https://doi.org/10.1093%2Fgenetics%2Fiyaa015)
blat <- read.table("raw_supp/FileS12_best_blatscore_xy_hit_length_filtered_uk_gwas.tsv", 
  sep="\t", header = TRUE, check.names=FALSE)
blat$per_identity[is.na(blat$per_identity)] <- 0 
snps_remaining <- blat %>% filter(per_identity < 90) #90% identity, >= 40bp overlap

#20,528 SNPs excluded, write to file
write.table(blat, "Kasimatis_removed_snps_BLAT.txt", 
  row.names = FALSE, sep="\t")


##### Additional site filtering (updated 5/16/24)
## Follows methods outlined in Ruzicka et al 2022 (remove sites with excess heterozygosity (FIS)
## sites with a deficit of minor allele homozygotes and significant missingness between males and females)

library(data.table)
library(dplyr)

#FIS, binomial, chisq functions
fis <- function(p,het){ (het/(2*p*(1-p))) - 1 }

binomial_pval <- function(x,n,p) {
  binom.test(x,n,p, alternative = c("less"))$p.value
}

chisq_missing <- function(x) {
  chisq.test(matrix(x,nrow = 2))$p.value
}

#load SNP and genotype data
males_geno <- fread("allchroms.males.gcount", sep="\t")
females_geno <- fread("allchroms.females.gcount", sep="\t")


#CGenotype counts

#Total genotype countd
males_geno$m_TOTAL <- males_geno$HOM_REF_CT+males_geno$HET_REF_ALT_CTS+males_geno$TWO_ALT_GENO_CTS
females_geno$f_TOTAL <- females_geno$HOM_REF_CT+females_geno$HET_REF_ALT_CTS+females_geno$TWO_ALT_GENO_CTS

# Heterozygote counts
males_geno$m_HET <- males_geno$HET_REF_ALT_CTS
females_geno$f_HET <- females_geno$HET_REF_ALT_CTS

#Homozygote counts for the reference allele
males_geno$m_HOM_REF <- males_geno$HOM_REF_CT
females_geno$f_HOM_REF <- females_geno$HOM_REF_CT

#Homozygote counts for the alternate allele
males_geno$m_HOM_ALT <- males_geno$TWO_ALT_GENO_CTS
females_geno$f_HOM_ALT <- females_geno$TWO_ALT_GENO_CTS

#Alternate allele frequency
males_geno$m_ALT_frq <- 1-((males_geno$HOM_REF_CT + 0.5*males_geno$HET_REF_ALT_CTS)/males_geno$m_TOTAL)
females_geno$f_ALT_frq <- 1-((females_geno$HOM_REF_CT + 0.5*females_geno$HET_REF_ALT_CTS)/females_geno$f_TOTAL)

# Reference allele frequency
males_geno$m_REF_frq <- 1-males_geno$m_ALT_frq
females_geno$f_REF_frq <- 1-females_geno$f_ALT_frq

#Minor allele frequency
males_geno$m_MAF <- ifelse(males_geno$m_ALT_frq <= 0.5, males_geno$m_ALT_frq, males_geno$m_REF_frq)
females_geno$f_MAF <- ifelse(females_geno$f_ALT_frq <= 0.5, females_geno$f_ALT_frq, females_geno$f_REF_frq)

#calculate missing
males_geno$m_MISSING <- 139865-males_geno$m_TOTAL
females_geno$f_MISSING <- 163959-females_geno$f_TOTAL

#### Combine sexes (x=male, y=female)
genotype_counts <- merge(males_geno, females_geno, by="ID")
genotype_counts$TOTAL <- genotype_counts$m_TOTAL + genotype_counts$f_TOTAL
genotype_counts$all_HET <- genotype_counts$m_HET + genotype_counts$f_HET
genotype_counts$all_HOM_REF <- genotype_counts$m_HOM_REF + genotype_counts$f_HOM_REF
genotype_counts$all_HOM_ALT <- genotype_counts$m_HOM_ALT + genotype_counts$f_HOM_ALT
genotype_counts$all_ALT_frq <- ((genotype_counts$m_ALT_frq*genotype_counts$m_TOTAL)+(genotype_counts$f_ALT_frq*genotype_counts$f_TOTAL))/(genotype_counts$m_TOTAL+genotype_counts$f_TOTAL) 
genotype_counts$MAF <- ifelse(( ((genotype_counts$m_ALT_frq*genotype_counts$m_TOTAL)+(genotype_counts$f_ALT_frq*genotype_counts$f_TOTAL))/(genotype_counts$m_TOTAL+genotype_counts$f_TOTAL)  )>0.5, 
                              1-( ((genotype_counts$m_ALT_frq*genotype_counts$m_TOTAL)+(genotype_counts$f_ALT_frq*genotype_counts$f_TOTAL))/(genotype_counts$m_TOTAL+genotype_counts$f_TOTAL)  ), 
                              ( ((genotype_counts$m_ALT_frq*genotype_counts$m_TOTAL)+(genotype_counts$f_ALT_frq*genotype_counts$f_TOTAL))/(genotype_counts$m_TOTAL+genotype_counts$f_TOTAL)  ))



#### Remove sites with excess heterozygosity (following Ruzicka et al, 2022)

genotype_counts$FIS <- fis((genotype_counts$m_ALT_frq+genotype_counts$f_ALT_frq)/2,(genotype_counts$all_HET/genotype_counts$TOTAL))
genotype_counts$m_FIS <- fis(genotype_counts$m_ALT_frq,(genotype_counts$m_HET/genotype_counts$m_TOTAL))
genotype_counts$f_FIS <- fis(genotype_counts$f_ALT_frq,(genotype_counts$f_HET/genotype_counts$f_TOTAL))

s <- 0.2 # Ruzicka et al = unrealistically strong selection coefficient at 0.2

genotype_counts$p_FIS <- p.adjust(pnorm(genotype_counts$FIS, mean=(1/(2*genotype_counts$TOTAL)) + (genotype_counts$MAF*(1-genotype_counts$MAF)/4)*(s/(1 - genotype_counts$MAF*s))^2, sd=sqrt(1/(genotype_counts$TOTAL)),lower.tail = F), "BH")

genotype_counts$p_FIS_m <- p.adjust(pnorm(genotype_counts$m_FIS, mean=(1/(2*genotype_counts$m_TOTAL)) + (genotype_counts$m_MAF*(1-genotype_counts$m_MAF)/4)*(s/(1 - genotype_counts$m_MAF*s))^2, sd=sqrt(1/(genotype_counts$m_TOTAL)),lower.tail = F), "BH")

genotype_counts$p_FIS_f <- p.adjust(pnorm(genotype_counts$f_FIS, mean=(1/(2*genotype_counts$f_TOTAL)) + (genotype_counts$f_MAF*(1-genotype_counts$f_MAF)/4)*(s/(1 - genotype_counts$f_MAF*s))^2, sd=sqrt(1/(genotype_counts$f_TOTAL)),lower.tail = F), "BH")


#### Remove sites with deficit of minor allele homozygotes 

#both
genotype_counts$p_MAH[genotype_counts$all_ALT_frq < 0.5] <- mapply(binomial_pval, #alt allele homozygotes
                                                                   genotype_counts$all_HOM_ALT[genotype_counts$all_ALT_frq<0.5],
                                                                   genotype_counts$TOTAL[genotype_counts$all_ALT_frq<0.5],
                                                                   genotype_counts$MAF[genotype_counts$all_ALT_frq<0.5]^2)

genotype_counts$p_MAH[genotype_counts$all_ALT_frq>=0.5] <- mapply(binomial_pval, #ref allele homozygotes
                                                                  genotype_counts$all_HOM_REF[genotype_counts$all_ALT_frq>=0.5],
                                                                  genotype_counts$TOTAL[genotype_counts$all_ALT_frq>=0.5],
                                                                  genotype_counts$MAF[genotype_counts$all_ALT_frq>=0.5]^2)
genotype_counts$p_MAH <- p.adjust(genotype_counts$p_MAH, "BH")


#males
genotype_counts$p_MAH_m[genotype_counts$m_ALT_frq < 0.5] <- mapply(binomial_pval, #alt allele homozygotes
                                                                   genotype_counts$m_HOM_ALT[genotype_counts$m_ALT_frq<0.5],
                                                                   genotype_counts$m_TOTAL[genotype_counts$m_ALT_frq<0.5],
                                                                   genotype_counts$m_MAF[genotype_counts$m_ALT_frq<0.5]^2)

genotype_counts$p_MAH_m[genotype_counts$m_ALT_frq>=0.5] <- mapply(binomial_pval, #ref allele homozygotes
                                                                  genotype_counts$m_HOM_REF[genotype_counts$m_ALT_frq>=0.5],
                                                                  genotype_counts$m_TOTAL[genotype_counts$m_ALT_frq>=0.5],
                                                                  genotype_counts$m_MAF[genotype_counts$m_ALT_frq>=0.5]^2)
genotype_counts$p_MAH_m <- p.adjust(genotype_counts$p_MAH_m, "BH")


#females
genotype_counts$p_MAH_f[genotype_counts$f_ALT_frq < 0.5] <- mapply(binomial_pval, #alt allele homozygotes
                                                                   genotype_counts$f_HOM_ALT[genotype_counts$f_ALT_frq<0.5],
                                                                   genotype_counts$f_TOTAL[genotype_counts$f_ALT_frq<0.5],
                                                                   genotype_counts$f_MAF[genotype_counts$f_ALT_frq<0.5]^2)

genotype_counts$p_MAH_f[genotype_counts$f_ALT_frq>=0.5] <- mapply(binomial_pval, #ref allele homozygotes
                                                                  genotype_counts$f_HOM_REF[genotype_counts$f_ALT_frq>=0.5],
                                                                  genotype_counts$f_TOTAL[genotype_counts$f_ALT_frq>=0.5],
                                                                  genotype_counts$f_MAF[genotype_counts$f_ALT_frq>=0.5]^2)
genotype_counts$p_MAH_f <- p.adjust(genotype_counts$p_MAH_f, "BH")


#### Filter for difference in missingness between males and females
genotype_counts$missing_pval <- apply(genotype_counts[,c("m_MISSING","m_TOTAL","f_MISSING","f_TOTAL")],1,function(x) chisq.test(matrix(x,nrow = 2))$p.value)
genotype_counts$missing_pval_bh <- p.adjust(genotype_counts$missing_pval, "fdr")

# Filter and export
genotype_counts_remove <- subset(genotype_counts, p_FIS < 0.05 | p_FIS_m < 0.05 | p_FIS_f < 0.05 | 
                                   p_MAH < 0.05 | p_MAH_m < 0.05 | p_MAH_f < 0.05 | missing_pval_bh < 0.05)

write.table(genotype_counts_remove$ID, 
            "genotype_filters_snps_toremove_5-23.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)