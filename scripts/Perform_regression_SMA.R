# Perform SMA regression with weightings on 27 traits
## 5/2/24

# libraries
library(data.table)
library(dplyr)

#directory
dr <- "/scratch/ukb/data/"    

#weighted RMA regression function 
sma_weighted <- function(x, y, sx, sy) {

  #weights
  wx <- sx^2
  wy <- sy^2
  v <- 1/(wx + wy)
  w <- v/sum(v)
    
  #weighted means  
  x_bar <- sum(w * x)
  y_bar <- sum(w * y)
    
  #weighted covariance and variances
  cov_xy <- sum(w * (x - x_bar) * (y - y_bar))
  var_x <- sum(w * (x - x_bar)^2)
  var_y <- sum(w * (y - y_bar)^2)
    
  #slope and intercept
  b1 <- sign(cov_xy)*(var_y / var_x)^0.5
  b0 <- y_bar - b1 * x_bar
    
  #SE
  n <- length(x)
  r <- cov_xy/sqrt(var_x*var_y)
  se<- abs(b1) * sqrt((1-r^2)/n)
    
  return(c(b0, b1,se))
}

#inputs
trait <- commandArgs(trailingOnly = TRUE)[1]
pexp <- commandArgs(trailingOnly = TRUE)[2]
pval_thres <- commandArgs(trailingOnly = TRUE)[3]

#trait file
filename_l <- paste0(dr, "traits/",trait,".gwas_ML_lowlfsr.txt")
GWAS_data_lfsr <- read.table(filename_l, header = TRUE)

#P-value threshold for GWAS snps
pval <- as.numeric(1)
#pval <- as.numeric(pval_thres)
GWAS_data_filtered_lfsr <- GWAS_data_lfsr[GWAS_data_lfsr$Male_lfsr <= pval | GWAS_data_lfsr$Female_lfsr <= pval,] 

#Effect size sex differences 
GWAS_data_filtered_lfsr$Effect_Diff <- (GWAS_data_filtered_lfsr$Female_pm - GWAS_data_filtered_lfsr$Male_pm)

#Effect size SEs
GWAS_data_filtered_lfsr$varEffect <- ((GWAS_data_filtered_lfsr$Female_psd)^2 + (GWAS_data_filtered_lfsr$Male_psd)^2)
GWAS_data_filtered_lfsr$sdEffect <- sqrt(GWAS_data_filtered_lfsr$varEffect)

#output summary GWAS file
write.table(GWAS_data_filtered_lfsr, 
  file = paste0(dr, "traits/",trait,".gwas_mash_summary.txt"), 
            sep = "\t",row.names = FALSE,append = FALSE)

################################################################################

#Bootstrap regression (viability)
slopes_V_l <- vector(length = 1000)

for(i in 1:1000) {
  #randomly sample 1 snp per block
  snp_sample_l <- GWAS_data_filtered_lfsr %>% group_by(Chrom,block) %>% slice_sample(n=1)
  
  #fit model (SMA and OLS)
  sma_via_l <- sma_weighted(snp_sample_l$Effect_Diff, snp_sample_l$s_v, 
    snp_sample_l$sdEffect, snp_sample_l$SE_v)
  slopes_V_l[i] <- sma_via_l[2]
}

# viability out
Z_V_l <- mean(slopes_V_l)/sd(slopes_V_l)


boot_via_l <- data.frame("Trait" = trait,
                         "Mode" = "Viability",
                         "Data" = "Lowest LFSR",
                         "mean_slope" = mean(slopes_V_l),
                         "SD_slope" = sd(slopes_V_l),
                         "Z" = Z_V_l,
                         "sample_size_regression" = nrow(snp_sample_l),
                         "SNPs_after_pval_filtering" = nrow(GWAS_data_filtered_lfsr))

# get output file and write it
write.table(boot_via_l, file = paste0(dr, "traits/",trait,".regression.lowestlfsr.mash.",pexp,".v.result"), 
  sep = "\t", row.names = FALSE,append = FALSE)


###########################################################################

#Bootstrap regression (fecundity)
slopes_F_l <- vector(length = 1000)

for(i in 1:1000) {
  #randomly sample 1 snp per block
  snp_sample_l <- GWAS_data_filtered_lfsr %>% group_by(Chrom,block) %>% slice_sample(n=1)
  
  #fit model (SMA and OLS)
  sma_fec_l <- sma_weighted(snp_sample_l$Effect_Diff, snp_sample_l$s_f, 
    snp_sample_l$sdEffect, snp_sample_l$SE_f)
  slopes_F_l[i] <- sma_fec_l[2]
}

# viability out
Z_F_l <- mean(slopes_F_l)/sd(slopes_F_l)


boot_fec_l <- data.frame("Trait" = trait,
                         "Mode" = "Fecundity",
                         "Data" = "Lowest LFSR",
                         "mean_slope" = mean(slopes_F_l),
                         "SD_slope" = sd(slopes_F_l),
                         "Z" = Z_F_l,
                         "sample_size_regression" = nrow(snp_sample_l),
                         "SNPs_after_pval_filtering" = nrow(GWAS_data_filtered_lfsr))

# get output file and write it
write.table(boot_fec_l, file = paste0(dr, "traits/",trait,".regression.lowestlfsr.mash.",pexp,".f.result"), 
  sep = "\t",row.names = FALSE,append = FALSE)


###########################################################################

#Bootstrap regression (total)
slopes_T_l <- vector(length = 1000)

for(i in 1:1000) {
  #randomly sample 1 snp per block
  snp_sample_l <- GWAS_data_filtered_lfsr %>% group_by(Chrom,block) %>% slice_sample(n=1)
  
  #fit model (SMA and OLS)
  sma_tot_l <- sma_weighted(snp_sample_l$Effect_Diff, snp_sample_l$s_t, snp_sample_l$sdEffect, snp_sample_l$SE_t)
  slopes_T_l[i] <- sma_tot_l[2]
}

# viability out
Z_T_l <- mean(slopes_T_l)/sd(slopes_T_l)


boot_tot_l <- data.frame("Trait" = trait,
                         "Mode" = "Total",
                         "Data" = "Lowest LFSR",
                         "mean_slope" = mean(slopes_T_l),
                         "SD_slope" = sd(slopes_T_l),
                         "Z" = Z_T_l,
                         "sample_size_regression" = nrow(snp_sample_l),
                         "SNPs_after_pval_filtering" = nrow(GWAS_data_filtered_lfsr))

# get output file and write it
write.table(boot_tot_l, file = paste0(dr, "traits/",trait,".regression.lowestlfsr.mash.",pexp,".t.result"), 
  sep = "\t",row.names = FALSE,append = FALSE)