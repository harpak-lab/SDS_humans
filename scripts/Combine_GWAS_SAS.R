#  Combine coefficients with MASH results
# updated - 3/15/24

# libraries
library(data.table)
library(dplyr)

#directory
dr <- "/scratch/ukb/data/"    

#current trait
trait <- commandArgs(trailingOnly = TRUE)[1]
sink(paste0("merge.",trait,".log"), type=c("output", "message"))
print("------------------------")
print(paste("Trait: ",trait))

# Load ML data (FILTERED)
file_ukb <- paste0(dr, "UKB_filtered_5-22.txt")
UKB_ML <- read.table(file_ukb, 
    sep="", header = TRUE, check.names=FALSE, stringsAsFactors = FALSE)

print(paste("Total windows: ", nrow(UKB_ML)))

# Associate leading snps in UKB ML dataset with haplotype blocks
blks <- paste0(dr, "UKB_snps_blocks.txt")
hapblock_snps <- read.table(blks, sep="", 
    header = TRUE, check.names=FALSE, stringsAsFactors = FALSE)
names(hapblock_snps)<- c("Chrom","Site1_Position","Site1_SNP","block")

#merge with data
UKB_ML <- merge(UKB_ML, hapblock_snps, by="Site1_SNP")

# Load GWAS data
filename_pval <- paste0(dr, "traits/", trait,".pvals.txt")
filename_effect <- paste0(dr, "traits/", trait,".effects.txt")

pval <- read.table(filename_pval, sep="", header = TRUE, check.names=FALSE)
effects <- read.table(filename_effect, sep="", header = TRUE, check.names=FALSE)

names(pval) <- c("Chrom","Position","MaleP","FemaleP")
names(effects) <- c("Chrom","Position","MaleB","MaleSE","FemaleB","FemaleSE")

#Load mash data (Zhu et al, 2023: https://doi.org/10.1016/j.xgen.2023.100297)
print("Loading mash data")
filename_mash_pm<- paste0(dr, "traits/", trait,"_mash_pm.txt")
filename_mash_psd<- paste0(dr, "traits/", trait, "_mash_psd.txt")
filename_mash_lfsr<- paste0(dr, "traits/", trait,"_mash_lfsr.txt")
mash_pm <-read.table(filename_mash_pm, sep="", header = TRUE, 
    check.names=FALSE)
mash_psd <-read.table(filename_mash_psd, sep="", header = TRUE, 
    check.names=FALSE)
mash_lfsr <- read.table(filename_mash_lfsr, sep="", header = TRUE, 
    check.names=FALSE)
names(mash_pm) <- c("Male_pm","Female_pm")
names(mash_psd) <- c("Male_psd","Female_psd")
names(mash_lfsr) <-  c("Male_lfsr","Female_lfsr")

#combine mash data into one
mash <- cbind(cbind(mash_pm,mash_psd),mash_lfsr)

#Merge mash with effects
effects <- cbind(effects,mash)

#Merge effect sizes with p-values
gwas_total <- merge(effects, pval,
 by=c("Chrom","Position"))

#Load in allele effects data (Zhu et al, 2023: https://doi.org/10.1016/j.xgen.2023.100297)
filename_male_effect <- paste0(dr, "traits/", trait,".male.effect.txt")
filename_female_effect <- paste0(dr, "traits/", trait,".female.effect.txt")
effect_males <-read.table(filename_male_effect, sep="", 
    header = FALSE, check.names=FALSE)
effect_females <- read.table(filename_female_effect, sep="", 
    header = FALSE, check.names=FALSE)
names(effect_males) <- c("Chrom","Position","SNP","A_male")
names(effect_females) <- c("Chrom","Position","SNP","A_female")

#merge into one effect file
allele_effects <- merge(effect_males,effect_females, 
    by=c("Chrom","Position","SNP"))

#Merge allele effect info with GWAS data
gwas_total_2 <- merge(gwas_total, allele_effects, 
    by=c("Chrom","Position"))

#Get unified dataframe 
UKB_ML_1 <- UKB_ML %>%
  left_join(gwas_total_2 %>% rename(Site1_SNP = SNP, 
                                    Site1_Position = Position,
                                    Site1_MaleB = MaleB,
                                    Site1_MaleSE = MaleSE,
                                    Site1_FemaleB = FemaleB,
                                    Site1_FemaleSE = FemaleSE,
                                    Site1_Male_pm = Male_pm,
                                    Site1_Female_pm = Female_pm,
                                    Site1_Male_psd = Male_psd,
                                    Site1_Female_psd = Female_psd,
                                    Site1_Male_lfsr = Male_lfsr,
                                    Site1_Female_lfsr = Female_lfsr,
                                    Site1_MaleP = MaleP,
                                    Site1_FemaleP = FemaleP,
                                    Site1_A_male = A_male,
                                    Site1_A_female = A_female),
            by = "Site1_SNP")

UKB_ML_2 <- UKB_ML_1 %>%
  left_join(gwas_total_2 %>% rename(Site0_SNP = SNP, 
                                    Site0_Position = Position,
                                    Site0_MaleB = MaleB,
                                    Site0_MaleSE = MaleSE,
                                    Site0_FemaleB = FemaleB,
                                    Site0_FemaleSE = FemaleSE,
                                    Site0_Male_pm = Male_pm,
                                    Site0_Female_pm = Female_pm,
                                    Site0_Male_psd = Male_psd,
                                    Site0_Female_psd = Female_psd,
                                    Site0_Male_lfsr = Male_lfsr,
                                    Site0_Female_lfsr = Female_lfsr,
                                    Site0_MaleP = MaleP,
                                    Site0_FemaleP = FemaleP,
                                    Site0_A_male = A_male,
                                    Site0_A_female = A_female),
            by = "Site0_SNP")

UKB_ML_2 <- UKB_ML_2[!(is.na(UKB_ML_2$Site1_Male_lfsr) & is.na(UKB_ML_2$Site0_Male_lfsr)), ]

UKB_ML_3 <- UKB_ML_2 %>% select(Chrom = Chrom.x, Window, 
                                Site1_Position = Site1_Position.x, Site1_SNP, Site1_MAF, 
                                Site1_Allele0, Site1_Allele1, Site1_MaleB,
                                Site1_MaleSE,Site1_FemaleB,Site1_FemaleSE,Site1_Male_pm,
                                Site1_Female_pm,Site1_Male_psd,Site1_Female_psd,
                                Site1_Male_lfsr,Site1_Female_lfsr,Site1_MaleP,Site1_FemaleP,
                                Site1_A_male,Site1_A_female,
                                Site0_Position = Site0_Position.x, Site0_SNP, Site0_MAF, 
                                Site0_Allele0, Site0_Allele1, Site0_MaleB,
                                Site0_MaleSE,Site0_FemaleB,Site0_FemaleSE,Site0_Male_pm,
                                Site0_Female_pm,Site0_Male_psd,Site0_Female_psd,
                                Site0_Male_lfsr,Site0_Female_lfsr,Site0_MaleP,Site0_FemaleP,
                                Site0_A_male,Site0_A_female,
                                s_v,SE_v,Z_v,s_f,SE_f,Z_f,s_t,SE_t,Z_t,
                                block, Site1_A_male, Site1_A_female, Site0_A_male, 
                                Site0_A_female)

### PICK SNP, LOWEST LFSR (viability)
print("Selecting Lowest LFSR SNPs.......")
large_value <- max(c(UKB_ML_3$Site1_Male_lfsr, UKB_ML_3$Site1_Female_lfsr, 
    UKB_ML_3$Site0_Male_lfsr, UKB_ML_3$Site0_Female_lfsr), na.rm = TRUE) + 99
lfsr_values <- cbind(UKB_ML_3$Site1_Male_lfsr, UKB_ML_3$Site1_Female_lfsr, 
    UKB_ML_3$Site0_Male_lfsr, UKB_ML_3$Site0_Female_lfsr)
lfsr_values[is.na(lfsr_values)] <- large_value

min_indices <- max.col(-lfsr_values, ties.method = "first")

SNP <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_SNP, UKB_ML_3$Site0_SNP)
Position <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Position, UKB_ML_3$Site0_Position)
SNP_allele_0 <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Allele0, UKB_ML_3$Site0_Allele0)
SNP_allele_1 <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Allele1, UKB_ML_3$Site0_Allele1)
Male_pm <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Male_pm, UKB_ML_3$Site0_Male_pm)
Male_psd <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Male_psd, UKB_ML_3$Site0_Male_psd)
Female_pm <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Female_pm, UKB_ML_3$Site0_Female_pm)
Female_psd <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Female_psd, UKB_ML_3$Site0_Female_psd)
A_male <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_A_male, UKB_ML_3$Site0_A_male)
A_female <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_A_female, UKB_ML_3$Site0_A_female)
Male_lfsr <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Male_lfsr, UKB_ML_3$Site0_Male_lfsr)
Female_lfsr <- ifelse(min_indices %in% c(1, 2), UKB_ML_3$Site1_Female_lfsr, UKB_ML_3$Site0_Female_lfsr)

# output
low_lfsr_snps <- data.frame("Chrom" = UKB_ML_3$Chrom, 
                            "Position" = Position,
                            "SNP" = SNP,
                            "SNP_allele_0" = SNP_allele_0,
                            "SNP_allele_1" = SNP_allele_1,
                            "s_v" = UKB_ML_3$s_v,
                            "SE_v" = UKB_ML_3$SE_v,
                            "s_f" = UKB_ML_3$s_f,
                            "SE_f" = UKB_ML_3$SE_f,
                            "s_t" = UKB_ML_3$s_t,
                            "SE_t" = UKB_ML_3$SE_t,
                            "Male_pm" = Male_pm,
                            "Male_psd" = Male_psd,
                            "Male_lfsr" = Male_lfsr,
                            "Female_pm" = Female_pm,
                            "Female_psd" = Female_psd,
                            "Female_lfsr" = Female_lfsr,
                            "A_male" = A_male,
                            "A_female" = A_female,
                            "block" = UKB_ML_3$block)

# Polarize data to be consistent with effect allele 
print("Polarizing......")
low_lfsr_snps$Male_pm <- ifelse(low_lfsr_snps$SNP_allele_1 == low_lfsr_snps$A_male, 
    ow_lfsr_snps$Male_pm, -(low_lfsr_snps$Male_pm))
low_lfsr_snps$Female_pm <- ifelse(low_lfsr_snps$SNP_allele_1 == low_lfsr_snps$A_female, 
    low_lfsr_snps$Female_pm, -(low_lfsr_snps$Female_pm))

#print log information
print(paste("SNPs in GWAS for trait): ", nrow(gwas_total_2)))
print(paste("Overlapping SNP windows in ML data (low lfsr): ", nrow(low_lfsr_snps)))
print(paste("Windows thrown out (low lfsr): ", (nrow(UKB_ML)-nrow(low_lfsr_snps))))

#Write out
print("------------------------")

write.table(low_lfsr_snps, file = paste0(trait,".gwas_ML_lowlfsr.txt"), 
    sep = "\t",row.names = FALSE,append = FALSE)
