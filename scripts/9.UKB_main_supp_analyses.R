############# Cole et al, 2024 Analyses (Figs 1C-3)
############# updated 4/25/24
##################################################
options(scipen = 999)

library(data.table) 
library(qqman)
library(ggplot2)
library(dplyr)

dr <- "/scratch/ukb/data/"

###################################################
## Simulation results

### Functions to calculate summary statistics

# MSE w/ bootstrapping
bootstrap_mse <- function(s_values, s, b) {
  B<-length(s_values)
  bootstrap_mses <- b
  
  for (b in 1:B) {
    bootstrap_sample <- sample(s_values, replace=TRUE)
    bootstrap_mses[b] <- mean((bootstrap_sample - s)^2)
  }
  
  lower_ci <- quantile(bootstrap_mses, probs = 0.025)
  upper_ci <- quantile(bootstrap_mses, probs = 0.975)
  
  cis <- c(lower_ci, upper_ci)
  
  return(cis)
}

# Calculate MSE and relative bias
getBias_MSE <- function(sims, st, px_hat, mse_boots) {
  
  sim_results <- do.call(rbind, lapply(sims, data.frame))
  sim_results$px_hat <- px_hat
  
  biases <- lapply(sims, function(sim) {
    list(
      Bias_int = mean(sim$int_s - st),
      Bias_sing = mean(sim$sing_s - st),
      abs_Bias_int = mean(abs(sim$int_s - st)),
      abs_Bias_sing = mean(abs(sim$sing_s - st)),
      RelErr_int = mean((sim$int_s / st) - 1),
      RelErr_sing = mean((sim$sing_s / st) - 1),
      abs_RelErr_int = mean(abs((sim$int_s / st) - 1)),
      abs_RelErr_sing = mean(abs((sim$sing_s / st) - 1))
    )
  })
  
  bias_err_values <- data.frame(
    d = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
    true_s = st,
    do.call(rbind, biases)
  )
  
  mse <- lapply(seq_along(sims), function(i) {
    sim_summary <- data.frame(sims[[i]])
    int_mse <- mean((sim_summary$int_s - st)^2)
    sing_mse <- mean((sim_summary$sing_s - st)^2)
    int_ci <- bootstrap_mse(sim_summary$int_s, st, mse_boots)
    sing_ci <- bootstrap_mse(sim_summary$sing_s, st, mse_boots)
    list(
      MSE_int = int_mse,
      MSE_sing = sing_mse,
      UpperCI_int = int_ci[[1]],
      LowerCI_int = int_ci[[2]],
      UpperCI_sing = sing_ci[[1]],
      LowerCI_sing = sing_ci[[2]]
    )
  })
  
  MSE <- data.frame(
    d = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
    true_s = st,
    do.call(rbind, mse)
  )
  
  res <- list(
    Results = sim_results,
    Bias = bias_err_values,
    MSE = MSE
  )
  
  return(res)
}


#####################
#Load sim data for three s values

sim_path_s1_p13 <- paste0(dr, "simulations/sim_d", 
                     0:5, 
                     "_s1_p13.sim")
sim_path_s01_p13 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s01_p13.sim")
sim_path_s001_p13 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s001_p13.sim")

sims_s1_p13 <- lapply(sim_path_s1_p13, fread, sep="\t", 
                      header=TRUE, check.names=FALSE)
sims_s01_p13 <- lapply(sim_path_s1_p13, fread, sep="\t", 
                       header=TRUE, check.names=FALSE)
sims_s001_p13 <- lapply(sim_path_s1_p13, fread, sep="\t", 
                        header=TRUE, check.names=FALSE)


#results tables (output all tables, bias, MSE calculations)
sim_s1_p13_result <- getBias_MSE(sims_s1_p13,st = 0.1,
                                 px_hat = 0.13, 
                                 mse_boots = 1000)
sim_s01_p13_result <- getBias_MSE(sims_s01_p13,st = 0.01,
                                  px_hat = 0.13, 
                                  mse_boots = 1000)
sim_s001_p13_result <- getBias_MSE(sims_s001_p13,st = 0.001,
                                   px_hat = 0.13, 
                                   mse_boots = 1000)

#export results tables

results_d<- list(sim_s1_p13_result$Results, 
                 sim_s01_p13_result$Results,
                 sim_s001_p13_result$Results)
bias_d<- list(sim_s1_p13_result$Bias, 
              sim_s01_p13_result$Bias,
                 sim_s001_p13_result$Bias)
mse_d<- list(sim_s1_p13_result$MSE, 
             sim_s01_p13_result$MSE,
                 sim_s001_p13_result$MSE)

#write out
names <- c("s1", "s01", "s001")
for (i in seq_along(results_d)) {
  write.table(results_d[[i]],
              file = paste0(dr, "simulations/sim_", 
                            names[i], "_p13_result.txt"),
              sep = "\t", row.names = FALSE)
}
for (i in seq_along(bias_d)) {
  write.table(bias_d[[i]],
              file = paste0(dr, "simulations/sim_", 
                            names[i], "_p13_bias.txt"),
              sep = "\t", row.names = FALSE)
}
for (i in seq_along(mse_d)) {
  write.table(mse_d[[i]],
              file = paste0(dr, "simulations/sim_", 
                            names[i], "_p13_MSE.txt"),
              sep = "\t", row.names = FALSE)
}

# Window filtering steps for UK Biobank
###################################################

#counts data
UKB_counts <- fread(paste0(dr, "all_Chroms.counts"))

#check all haps present
missing_haps <- UKB_counts %>% 
  group_by(Leading_SNP) %>%
  summarize(MissingHaplotypes = setdiff(c(0, 1, 2, 3), unique(haplotype))) %>%
  filter(length(MissingHaplotypes) > 0) %>%
  pull(Leading_SNP)

##extract low hap counts 
filter_UKB <- UKB_counts %>% 
  filter(female_counts < 10 & male_counts < 10)

# Return a vector of Leading_SNPs from the filtered dataframe
snp_list_remove <- c(unique(filter_UKB$Leading_SNP),unique(missing_haps))


#empirical results
UKB_via <- fread(paste0(dr, "allchroms2.v.result"))
UKB_fec <- fread(paste0(dr, "allchroms2.f.result")) 
UKB_tot <- fread(paste0(dr, "allchroms2.t.result"))

names(UKB_tot)[13]<-"r_10_t"
names(UKB_tot)[16]<-"D_10_t"
names(UKB_tot)[31]<-"Outside_boot_tot"

#permuted results
UKB_via_perm <- fread(paste0(dr, "allchroms2.permuted.v.result"))
UKB_tot_perm <- fread(paste0(dr, "allchroms2.permuted.t.result"))
UKB_fec_perm <- fread(paste0(dr, "allchroms2.permuted.f.result"))

names(UKB_tot_perm)[13]<-"r_10_t"
names(UKB_tot_perm)[16]<-"D_10_t"
names(UKB_tot_perm)[31]<-"Outside_boot_tot"

#merge all data
UKB_obs <- merge(merge(UKB_via, UKB_fec, by = intersect(names(UKB_via), names(UKB_fec))),
                 UKB_tot, by = intersect(names(UKB_via), names(UKB_tot)))
UKB_perm <- merge(merge(UKB_via_perm, UKB_fec_perm, by = intersect(names(UKB_via_perm), names(UKB_fec_perm))),
                 UKB_tot_perm, by = intersect(names(UKB_via_perm), names(UKB_tot_perm))) %>% 
  select(Chrom, Window, Site1_Position, Site1_SNP, Site1_MAF ,Site1_Allele0, #set column names for permuted data
         Site1_Allele1,Site0_Position, Site0_SNP, 
         Site0_MAF ,Site0_Allele0, Site0_Allele1,
         r_10_perm = r_10, r_1x_perm = r_1x, r_x0_perm = r_x0, 
         D_10_perm = D_10, D_1x_perm = D_1x, D_x0_perm = D_x0,
         s_v_perm = s_int_viability, ML_s_v_perm=ML_int_viability, 
         SE_v_perm = SE_v,Z_v_perm = Z_v,s_v_pval_LRT_perm = Pval_int_lrt_viability,
         s_v_pval_perm=Pval_int_boot_viability, 
         s_site1_viability_perm = s_site1_viability, ML_site1_s_v_perm = ML_site1_viability, 
         Pval_site1_lrt_viability_perm = Pval_site1_lrt_viability,s_site0_viability_perm=s_site0_viability,
         ML_site0_s_v_perm = ML_site0_viability,Pval_site0_lrt_viability_perm = Pval_site0_lrt_viability, 
         Outside_boot_via_perm = Outside_boot_via,
         s_f_perm = s_int_fecundity,SE_f_perm = SE_f, Z_f_perm = Z_f, 
         s_f_pval_perm=Pval_int_boot_fecundity, s_site1_fecundity_perm = s_site1_fecundity, 
         s_site0_fecundity_perm=s_site0_fecundity,
         r_10_t_perm = r_10_t, r_1x_t_perm = r_1x_t, r_x0_t_perm = r_x0_t, 
         D_10_t_perm = D_10_t, D_1x_t_perm = D_1x_t, D_x0_t_perm = D_x0_t,
         s_t_perm = s_int_total, ML_s_t_perm=ML_int_total, 
         SE_t_perm = SE_t, Z_t_perm = Z_t, s_t_pval_LRT_perm = Pval_int_lrt_total,
         s_t_pval_perm=Pval_int_boot_total, s_site1_total_perm = s_site1_total, 
         ML_site1_s_t_perm = ML_site1_total, 
         Pval_site1_lrt_total_perm = Pval_site1_lrt_total, 
         s_site0_total_perm= s_site0_total, ML_site0_s_t_perm = ML_site0_total,
         Pval_site0_lrt_total_perm = Pval_site0_lrt_total, 
         Outside_boot_tot_perm = Outside_boot_tot)

UKB_all <- merge(UKB_obs, UKB_perm, by = intersect(names(UKB_obs), names(UKB_perm)))

#filter out low counts
UKB_via_filt1 <- UKB_via[!UKB_via$Site1_SNP %in% snp_list_remove, ]
UKB_via_filt1 <- UKB_via_filt1 %>% filter(abs(r_10) < 1.0)

UKB_via_perm_filt1 <- UKB_via_perm[!UKB_via_perm$Site1_SNP %in% snp_list_remove, ]
UKB_via_perm_filt1 <- UKB_via_perm_filt1 %>% filter(abs(r_10) < 1.0)


## filter out r_10 >= 0.1
UKB_via_filt2 <- UKB_via_filt1 %>% filter(abs(r_10) >= 0.1)
UKB_via_perm_filt2 <- UKB_via_perm_filt1 %>% filter(abs(r_10) >= 0.1)


#remove MHC
UKB_via_filt2  <- UKB_via_filt2  %>%
  filter(!(Chrom == 6 & Site1_Position >= 28477797 & Site1_Position <= 33448354))
UKB_via_perm_filt2  <- UKB_via_perm_filt2  %>%
  filter(!(Chrom == 6 & Site1_Position >= 28477797 & Site1_Position <= 33448354))

#ensure same windows
UKB_via_perm_filt2 <- UKB_via_perm_filt2[UKB_via_perm_filt2$Site1_SNP %in% UKB_via_filt2$Site1_SNP, ]
UKB_via_filt2 <- UKB_via_filt2[UKB_via_filt2$Site1_SNP %in% UKB_via_perm_filt2$Site1_SNP, ]

## filter total and fecundity
UKB_fec_filt2 <- UKB_fec[UKB_fec$Site1_SNP %in% UKB_via_filt2$Site1_SNP, ]
UKB_tot_filt2 <- UKB_tot[UKB_tot$Site1_SNP %in% UKB_via_filt2$Site1_SNP, ]

UKB_fec_perm_filt2 <- UKB_fec_perm[UKB_fec_perm$Site1_SNP %in% UKB_via_perm_filt2$Site1_SNP, ]
UKB_tot_perm_filt2 <- UKB_tot_perm[UKB_tot_perm$Site1_SNP %in% UKB_via_perm_filt2$Site1_SNP, ]

### merge observed with permuted
UKB_via_merged <- merge(UKB_via_filt2,UKB_via_perm_filt2, by="Site1_SNP")

UKB_tot_merged <- merge(UKB_tot_filt2,UKB_tot_perm_filt2, by="Site1_SNP")

UKB_fec_merged <- merge(UKB_fec_filt2,UKB_fec_perm_filt2, by="Site1_SNP")


UKB_via_merged<- UKB_via_merged %>% select(Chrom = Chrom.x, Window = Window.x, 
                                           Site1_Position=Site1_Position.x, Site1_SNP, Site1_MAF=Site1_MAF.x,Site1_Allele0=Site1_Allele0.x, 
                                           Site1_Allele1=Site1_Allele1.x,Site0_Position=Site0_Position.x, Site0_SNP=Site0_SNP.x, 
                                           Site0_MAF=Site0_MAF.x,Site0_Allele0=Site0_Allele0.x, Site0_Allele1=Site0_Allele1.x,
                                           s_v = s_int_viability.x, ML_s_v=ML_int_viability.x, s_v_pval=Pval_int_boot_viability.x, 
                                           s_v_pval_LRT = Pval_int_lrt_viability.x, SE_v = SE_v.x, Z_v = Z_v.x,
                                           s_site1_viability= s_site1_viability.x, ML_site1_s_v = ML_site1_viability.x, 
                                           s_site0_viability=s_site0_viability.x, ML_site0_s_v = ML_site0_viability.x, 
                                           s_v_perm = s_int_viability.y, ML_s_v_perm=ML_int_viability.y, 
                                           s_v_pval_perm=Pval_int_boot_viability.y, s_v_pval_LRT_perm = Pval_int_lrt_viability.y, 
                                           SE_v_perm = SE_v.y, Z_v_perm = Z_v.y, s_site1_viability_perm = s_site1_viability.y, 
                                           ML_site1_s_v_perm = ML_site1_viability.y, s_site0_viability_perm=s_site0_viability.y, 
                                           ML_site0_s_v_perm = ML_site0_viability.y)

UKB_tot_merged<- UKB_tot_merged %>% select(Chrom = Chrom.x, Window = Window.x,
                                           Site1_Position=Site1_Position.x, Site1_SNP, Site1_MAF=Site1_MAF.x,Site1_Allele0=Site1_Allele0.x,
                                           Site1_Allele1=Site1_Allele1.x, Site0_Position=Site0_Position.x, Site0_SNP=Site0_SNP.x,
                                           Site0_MAF=Site0_MAF.x,Site0_Allele0=Site0_Allele0.x, Site0_Allele1=Site0_Allele1.x,
                                           s_t = s_int_total.x, ML_s_t=ML_int_total.x, s_t_pval=Pval_int_boot_total.x,
                                           s_t_pval_LRT = Pval_int_lrt_total.x, SE_t = SE_t.x, Z_t = Z_t.x,
                                           s_site1_total= s_site1_total.x, ML_site1_s_t = ML_site1_total.x,
                                           s_site0_total=s_site0_total.x, ML_site0_s_t = ML_site0_total.x,
                                           s_t_perm = s_int_total.y, ML_s_t_perm=ML_int_total.y, s_t_pval_perm=Pval_int_boot_total.y,
                                           s_t_pval_LRT_perm = Pval_int_lrt_total.y, SE_t_perm = SE_t.y, Z_t_perm = Z_t.y,
                                           s_site1_total_perm = s_site1_total.y, ML_site1_s_t_perm = ML_site1_total.y,
                                           s_site0_total_perm=s_site0_total.y, ML_site0_s_t_perm = ML_site0_total.y)


UKB_fec_merged<- UKB_fec_merged %>% select(Chrom = Chrom.x, Window = Window.x, 
                                           Site1_Position=Site1_Position.x, Site1_SNP, Site1_MAF=Site1_MAF.x,Site1_Allele0=Site1_Allele0.x, 
                                           Site1_Allele1=Site1_Allele1.x, Site0_Position=Site0_Position.x, Site0_SNP=Site0_SNP.x, 
                                           Site0_MAF=Site0_MAF.x,Site0_Allele0=Site0_Allele0.x, Site0_Allele1=Site0_Allele1.x,
                                           s_f = s_int_fecundity.x, s_f_pval=Pval_int_boot_fecundity.x, 
                                           SE_f = SE_f.x, Z_f = Z_f.x,
                                           s_site1_fecundity= s_site1_fecundity.x,  
                                           s_site0_fecundity=s_site0_fecundity.x, 
                                           s_f_perm = s_int_fecundity.y, s_f_pval_perm=Pval_int_boot_fecundity.y, 
                                           SE_f_perm = SE_f.y, Z_f_perm = Z_f.y,
                                           s_site1_fecundity_perm = s_site1_fecundity.y, 
                                           s_site0_fecundity_perm=s_site0_fecundity.y)


# Load LD pruned SNPs
UKB_pruned <- readLines(paste0(dr, "all_chrs_pruned_02.txt"))

UKB_via_merged_pruned <- UKB_via_merged[UKB_via_merged$Site1_SNP %in% UKB_pruned, ]
UKB_fec_merged_pruned <- UKB_fec_merged[UKB_fec_merged$Site1_SNP %in% UKB_pruned, ]
UKB_tot_merged_pruned <- UKB_tot_merged[UKB_tot_merged$Site1_SNP %in% UKB_pruned, ]

### Merge all data
UKB_filt <- merge(merge(UKB_via_merged, UKB_fec_merged, by="Site1_SNP"),UKB_tot_merged, by="Site1_SNP") %>%
  select(-contains(c('.x','.y'))) %>% select(Chrom, Window, Site1_Position, Site1_SNP, Site1_MAF,
                   Site1_Allele0,Site1_Allele1, Site0_Position, 
                   Site0_SNP,Site0_MAF,Site0_Allele0,Site0_Allele1,
                   everything())
UKB_filt <-UKB_filt[order(UKB_filt$Chrom, UKB_filt$Window),]


UKB_filt_pruned <- merge(merge(UKB_via_merged_pruned, UKB_fec_merged_pruned, by="Site1_SNP"),
                         UKB_tot_merged_pruned, by="Site1_SNP") %>%
  select(-contains(c('.x','.y'))) %>% select(Chrom, Window, Site1_Position, Site1_SNP, Site1_MAF,
                                             Site1_Allele0,Site1_Allele1, Site0_Position, 
                                             Site0_SNP,Site0_MAF,Site0_Allele0,Site0_Allele1,
                                             everything())
UKB_filt_pruned <-UKB_filt_pruned[order(UKB_filt_pruned$Chrom, UKB_filt_pruned$Window),]


#export
write.table(UKB_all, 
            file = paste0(dr, "UKB_all_5-22.txt"), 
            sep = "\t", row.names = FALSE)
write.table(UKB_filt, 
            file = paste0(dr, "UKB_filtered_5-22.txt"), 
            sep = "\t", row.names = FALSE)
write.table(UKB_filt_pruned, 
            file = paste0(dr, "UKB_filtered_pruned_5-22.txt"), 
            sep = "\t", row.names = FALSE)


write.table(UKB_all,"UKB_all_5-22.txt", 
            sep = "\t", row.names = FALSE)
write.table(UKB_filt,"UKB_filtered_5-22.txt", 
            sep = "\t", row.names = FALSE)
write.table(UKB_filt_pruned, "UKB_filtered_pruned_5-22.txt", 
            sep = "\t", row.names = FALSE)



#################** Statistical Test, Mann Whitney/KS, observed vs permuted

##### Viability
#no pruning
summary(abs(UKB_filt$s_v))
summary(abs(UKB_filt$s_v_perm))

wilcox.test(abs(UKB_filt$s_v), abs(UKB_filt$s_v_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt$s_v), abs(UKB_filt$s_v_perm))

#with pruning
summary(abs(UKB_filt_pruned$s_v))
summary(abs(UKB_filt_pruned$s_v_perm))

wilcox.test(abs(UKB_filt_pruned$s_v), abs(UKB_filt_pruned$s_v_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt_pruned$s_v), abs(UKB_filt_pruned$s_v_perm))

##### Fecundity
#no pruning
summary(abs(UKB_filt$s_f))
summary(abs(UKB_filt$s_f_perm))

wilcox.test(abs(UKB_filt$s_f), abs(UKB_filt$s_f_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt$s_f), abs(UKB_filt$s_f_perm))

#with pruning
summary(abs(UKB_filt_pruned$s_f))
summary(abs(UKB_filt_pruned$s_f_perm))
        
wilcox.test(abs(UKB_filt_pruned$s_f), abs(UKB_filt_pruned$s_f_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt_pruned$s_f), abs(UKB_filt_pruned$s_f_perm))

##### Total
#no pruning
summary(abs(UKB_filt$s_t))
summary(abs(UKB_filt$s_t_perm))

wilcox.test(abs(UKB_filt$s_t), abs(UKB_filt$s_t_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt$s_t), abs(UKB_filt$s_t_perm))

#with pruning
summary(abs(UKB_filt_pruned$s_t))
summary(abs(UKB_filt_pruned$s_t_perm))

wilcox.test(abs(UKB_filt_pruned$s_v), abs(UKB_filt_pruned$s_v_perm), 
            alternative = "greater")
ks.test(abs(UKB_filt_pruned$s_v), abs(UKB_filt_pruned$s_v_perm))


#################* Binned enrichment analysis

## chi-squared test
chi_sq <- function(r, index) {
  obs <- r$obs_count[index]
  perm <-r$perm_count[index]
  total_obs <- sum(r$obs_count) - obs
  total_perm<- sum(r$perm_count) - perm
  
  matrix <- matrix(c(obs, total_obs, perm, total_perm), 
                   nrow = 2)
  tst <- chisq.test(matrix)
  
  return(tst$p.value)
}


##### VIABILITY

#create separate dfs
s_v_observed <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_v)
s_v_permuted <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_v_perm)
s_v_observed$s_v_obs <- abs(s_v_observed$s_v)
s_v_permuted$s_v_perm <- abs(s_v_permuted$s_v_perm)

# binning
bin_edges <- quantile(s_v_permuted$s_v_perm, 
                      probs = seq(0, 1, length.out = 12), 
                      na.rm = TRUE)
s_v_permuted$bin <- cut(s_v_permuted$s_v_perm, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin medians
bin_medians <- s_v_permuted %>%
  group_by(bin) %>%
  summarise(median_val = median(s_v_perm, na.rm = TRUE)) %>%
  mutate(bin_label = as.character(median_val))  
bin_medians <- as.data.frame(bin_medians)

# labels
s_v_permuted <- s_v_permuted %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

# Bin
s_v_observed$bin <- cut(s_v_observed$s_v_obs, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin_label to s_v_observed
s_v_observed <- s_v_observed %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

obs_counts <- s_v_observed %>%
  group_by(bin) %>%
  summarise(obs_count = n())

obs_counts <- na.omit(obs_counts)

# counts
perm_counts <- s_v_permuted %>%
  group_by(bin) %>%
  summarise(perm_count = n())

# Merge
via_results_pruned<- merge(obs_counts, 
                           perm_counts, 
                           by = "bin", 
                           all = TRUE)
via_results_pruned$obs_null <- via_results_pruned$obs_count-via_results_pruned$perm_count
via_results_pruned$obs_null_percent <- (via_results_pruned$obs_null / 
                                          via_results_pruned$obs_count) * 100

# get chi squared p-values
via_results_pruned$p_value <- sapply(1:nrow(via_results_pruned), 
                                     function(x) chi_sq(via_results_pruned, x))
via_results_pruned$bin_label <- round(as.numeric(via_results_pruned$bin)*100, 
                                      digits=2)



##### FECUNDITY

# create separate dfs
s_f_observed <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_f)
s_f_permuted <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_f_perm)
s_f_observed$s_f_obs <- abs(s_f_observed$s_f)
s_f_permuted$s_f_perm <- abs(s_f_permuted$s_f_perm)

# binning
bin_edges <- quantile(s_f_permuted$s_f_perm, 
                      probs = seq(0, 1, length.out = 12), 
                      na.rm = TRUE)
s_f_permuted$bin <- cut(s_f_permuted$s_f_perm, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin medians
bin_medians <- s_f_permuted %>%
  group_by(bin) %>%
  summarise(median_val = median(s_f_perm, 
                                na.rm = TRUE)) %>%
  mutate(bin_label = as.character(median_val)) 

bin_medians <- as.data.frame(bin_medians)

# labels
s_f_permuted <- s_f_permuted %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

# bin
s_f_observed$bin <- cut(s_f_observed$s_f_obs, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin_label to s_f_observed 
s_f_observed <- s_f_observed %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

obs_counts <- s_f_observed %>%
  group_by(bin) %>%
  summarise(obs_count = n())
obs_counts <- na.omit(obs_counts)

# counts
perm_counts <- s_f_permuted %>%
  group_by(bin) %>%
  summarise(perm_count = n())

# Merge
fec_results_pruned<- merge(obs_counts, 
                           perm_counts, by = "bin", 
                           all = TRUE)
fec_results_pruned$obs_null <- fec_results_pruned$obs_count-fec_results_pruned$perm_count
fec_results_pruned$obs_null_percent <- (fec_results_pruned$obs_null / 
                                          fec_results_pruned$obs_count) * 100

# get chi squared p-values
fec_results_pruned$p_value <- sapply(1:nrow(fec_results_pruned), 
                                     function(x) chi_sq(fec_results_pruned, x))
fec_results_pruned$bin_label <- round(as.numeric(fec_results_pruned$bin)*100, 
                                      digits=2)


##### TOTAL

# create separate dfs
s_t_observed <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_t)
s_t_permuted <- UKB_filt_pruned %>% 
  select(Chrom,Site1_SNP,s_t_perm)
s_t_observed$s_t_obs <- abs(s_t_observed$s_t)
s_t_permuted$s_t_perm <- abs(s_t_permuted$s_t_perm)

# binning
bin_edges <- quantile(s_t_permuted$s_t_perm, 
                      probs = seq(0, 1, length.out = 12), 
                      na.rm = TRUE)
s_t_permuted$bin <- cut(s_t_permuted$s_t_perm, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin medians
bin_medians <- s_t_permuted %>%
  group_by(bin) %>%
  summarise(median_val = median(s_t_perm, na.rm = TRUE)) %>%
  mutate(bin_label = as.character(median_val))  
bin_medians <- as.data.frame(bin_medians)

# labels
s_t_permuted <- s_t_permuted %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

# bin
s_t_observed$bin <- cut(s_t_observed$s_t_obs, 
                        breaks = bin_edges, 
                        include.lowest = TRUE, 
                        labels = FALSE)

# bin_label to s_t_observed 
s_t_observed <- s_t_observed %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) 

obs_counts <- s_t_observed %>%
  group_by(bin) %>%
  summarise(obs_count = n())

obs_counts <- na.omit(obs_counts)

# counts
perm_counts <- s_t_permuted %>%
  group_by(bin) %>%
  summarise(perm_count = n())

# Merge 
tot_results_pruned<- merge(obs_counts, perm_counts, 
                           by = "bin", 
                           all = TRUE)
tot_results_pruned$obs_null <- tot_results_pruned$obs_count-tot_results_pruned$perm_count
tot_results_pruned$obs_null_percent <- (tot_results_pruned$obs_null / 
                                          tot_results_pruned$obs_count) * 100


# get chi squared p-values
tot_results_pruned$p_value <- sapply(1:nrow(tot_results_pruned), 
                                     function(x) chi_sq(tot_results_pruned, x))
tot_results_pruned$bin_label <- round(as.numeric(tot_results_pruned$bin)*100, 
                                      digits=2)


#export
write.table(via_results_pruned, 
            file = paste0(dr, "UKB_v_bins_5-20.txt"), 
            sep = "\t", row.names = FALSE)
write.table(fec_results_pruned, 
            file = paste0(dr, "UKB_f_bins_5-20.txt"), 
            sep = "\t", row.names = FALSE)
write.table(tot_results_pruned, 
            file = paste0(dr, "UKB_t_bins_5-20.txt"), 
            sep = "\t", row.names = FALSE)



#################* Significant windows

UKB_via_filt2$FDR_boot <- p.adjust(UKB_via_filt2$Pval_int_boot_viability, 
                                   method = "hochberg")
UKB_fec_filt2$FDR_boot <- p.adjust(UKB_fec_filt2$Pval_int_boot_fecundity, 
                                   method = "hochberg")
UKB_tot_filt2$FDR_boot <- p.adjust(UKB_tot_filt2$Pval_int_boot_total, 
                                   method = "hochberg")

UKB_via_sig_windows <- UKB_via_filt2[UKB_via_filt2$FDR_boot < 0.05, ]
UKB_fec_sig_windows <- UKB_fec_filt2[UKB_fec_filt2$FDR_boot < 0.05, ]
UKB_tot_sig_windows <- UKB_tot_filt2[UKB_tot_filt2$FDR_boot < 0.05, ]

sig_vars_via <- unique(c(UKB_via_sig_windows$Site1_SNP, 
                         UKB_via_sig_windows$Site0_SNP))
sig_vars_tot <- unique(c(UKB_tot_sig_windows$Site1_SNP,
                         UKB_tot_sig_windows$Site0_SNP))

write(sig_vars_via,file = paste0(dr, "UKB_via_FDR_wins_5-20.txt"))
write(sig_vars_via,file = paste0(dr, "UKB_tot_FDR_wins_5-20.txt"))


#################** Viability vs fecundity 

#correlation and lm (coefficients)
cor.test(UKB_filt_pruned$s_f,UKB_filt_pruned$s_v)

vf_lm_s_pruned <- lm(s_v~s_f, data=UKB_filt_pruned)
summary(vf_lm_s_pruned)


#correlation and lm (Z-scores)
cor.test(UKB_filt_pruned$Z_f,UKB_filt_pruned$Z_v)

vf_lm_Z_pruned <- lm(Z_v~Z_f, data=UKB_filt_pruned)
summary(vf_lm_Z_pruned)



##################** UK Biobank Allele Freqs/counts from hap counts

#function to get allele frequencies 
get_allele_freqs <- function(data, SNP){
  
  #output DF
  output <- data.frame()
  
  #get window data
  window_counts <- data[data$Leading_SNP == SNP, ]
  Chrom <- unique(window_counts$Chrom)
  Window <- unique(window_counts$Window)
  SNP1 <- unique(window_counts$Leading_SNP)
  SNP1_Position <- unique(window_counts$Leading_Position)
  SNP1_MAF <- unique(window_counts$Leading_MAF)
  
  #get counts
  NF <- sum(window_counts$female_counts)
  NM <- sum(window_counts$male_counts)
  site1_counts <- window_counts[window_counts$haplotype %in% c(2, 3),]
  site0_counts <- window_counts[window_counts$haplotype %in% c(1, 3),]
  nF1_site1 <-  sum(site1_counts$female_counts)
  nF0_site1 <-  NF - nF1_site1
  nM1_site1 <-  sum(site1_counts$male_counts)
  nM0_site1 <-  NM - nM1_site1
  nF1_site0 <-  sum(site0_counts$female_counts)
  nF0_site0 <-  NF - nF1_site0
  nM1_site0 <-  sum(site0_counts$male_counts)
  nM0_site0 <-  NM - nM1_site0
  nF1_f_site1 <-  sum(site1_counts$adj_female_counts)
  nF0_f_site1 <-  NF - nF1_f_site1 
  nM1_f_site1 <-  sum(site1_counts$adj_male_counts)
  nM0_f_site1 <-  NM - nM1_f_site1
  nF1_f_site0 <-  sum(site0_counts$adj_female_counts)
  nF0_f_site0 <-  NF - nF1_f_site0 
  nM1_f_site0 <-  sum(site0_counts$adj_male_counts)
  nM0_f_site0 <-  NM - nM1_f_site0
  
  #allele freqs
  p1 <- (nF1_site1+nM1_site1)/(NF+NM)
  p1_female <- nF1_site1/NF
  p1_male <- nM1_site1/NM
  p1_sexavg <- (p1_female+p1_male)/2
  p1_prime_f <- nF1_f_site1/NF
  p1_prime_m <- nM1_f_site1/NM
  p0 <- (nF1_site0+nM1_site0)/(NF+NM)
  p0_female <-nF1_site0/NF
  p0_male <- nM1_site0/NM
  p0_sexavg <- (p0_female+p0_male)/2
  p0_prime_f <- nF1_f_site0/NF
  p0_prime_m <- nM1_f_site0/NM
  
  #Output results 
  out <- data.frame(Chrom,Window,SNP1_Position,SNP1,SNP1_MAF,
                    nF1_site1, nF0_site1, nM1_site1, nM0_site1,
                    nF1_site0, nF0_site0, nM1_site0, nM0_site0,
                    nF1_f_site1, nF0_f_site1, nM1_f_site1, nM0_f_site1,
                    nF1_f_site0, nF0_f_site0, nM1_f_site0, nM0_f_site0,
                    p1,p1_female,p1_male,p1_sexavg,p1_prime_f,p1_prime_m,
                    p0,p0_female,p0_male,p0_sexavg,p0_prime_f,p0_prime_m)
  
  output <- rbind(output, out)
  
  return(output)
}

#obs and permuted counts functions
get_allele_freqs_obs <- function(SNPs) {
  SNP <- SNPs
  get_allele_freqs(data = hap_counts, 
                   SNP=SNP)
}

get_allele_freqs_perm <- function(SNPs) {
  SNP <- SNPs
  get_allele_freqs(data = hap_counts_perm, 
                   SNP=SNP)
}


#function to estimate per-site ML S_hats analytically

via_s <- function(nF1,nF0,nM1,nM0){  #viability
  NM <- nM0 + nM1
  NF <- nF0 + nF1
  a <- 2*NF*NM*((nF1*nM0)-(nF0*nM1))
  b <- (2*nF0*nF1*(NM)^2) + (2*nM0*nM1*(NF)^2) + ((nF0*nM1)-(nF1*nM0))^2
  s_v <- a/b
  return(s_v)
}

tot_s <- function(p1_f,p1_m,p1_avg){  #total
  a <- (p1_f-p1_m)
  b <- 2*p1_avg*(1-p1_avg)
  s_t <- a/b
  return(s_t)
}

fec_s <- function(sv,st,p){
  a <- (st-sv)
  b <- (1+(1-p)*sv)*(1-p*sv)
  s_f <- a/b
  return(s_f)
}


##load counts data
hap_counts <- UKB_counts
hap_counts_perm <- read.table("allChroms.permuted.counts", sep="", 
                              header = TRUE, 
                              check.names=FALSE, 
                              stringsAsFactors = FALSE)
hap_counts_perm <- hap_counts_perm[hap_counts_perm$Leading_SNP %in% hap_counts$Leading_SNP,]

#unique chrom+windows
cw <- unique(hap_counts[, "Leading_SNP"])
cw_list <- split(cw, seq(length(cw)))
cw_p <- unique(hap_counts_perm[, "Leading_SNP"])
cw_list_p <- split(cw_p, seq(length(cw_p)))

#get freqs
allele_freqs_obs <- lapply(cw_list, get_allele_freqs_obs)
allele_freqs_obs <- do.call(rbind, allele_freqs_obs)
allele_freqs_perm <- lapply(cw_list_p, get_allele_freqs_perm)
allele_freqs_perm <- do.call(rbind, allele_freqs_perm)

#x=obs, y=perm
UKB_afs<-merge(allele_freqs_obs,allele_freqs_perm, 
               by=c("Chrom","Window","SNP1_Position","SNP1","SNP1_MAF"))
rename_suff <- function(x) {
  names(x) <- sub("\\.x$", ".obs", names(x))
  names(x) <- sub("\\.y$", ".perm", names(x))
  return(x)
}
UKB_afs <- rename_suff(UKB_afs)

#export
write.table(UKB_afs, file = paste0(dr, "UKB_allele_frqs_5-20.txt"), 
            sep = "\t",
            row.names = FALSE,
            append = FALSE)


#calculate per-site selection coefficients analytically

UKB_afs[ , 6:ncol(UKB_afs)] <- lapply(UKB_afs[ , 6:ncol(UKB_afs)], as.numeric)

UKB_afs$s_site1_v_obs<- via_s(UKB_afs$nF1_site1.obs,UKB_afs$nF0_site1.obs,
                              UKB_afs$nM1_site1.obs,UKB_afs$nM0_site1.obs)
UKB_afs$s_site1_v_perm<- via_s(UKB_afs$nF1_site1.perm,UKB_afs$nF0_site1.perm,
                              UKB_afs$nM1_site1.perm,UKB_afs$nM0_site1.perm)
UKB_afs$s_site0_v_obs<- via_s(UKB_afs$nF1_site0.obs,UKB_afs$nF0_site0.obs,
                              UKB_afs$nM1_site0.obs,UKB_afs$nM0_site0.obs)
UKB_afs$s_site0_v_perm<- via_s(UKB_afs$nF1_site0.perm,UKB_afs$nF0_site0.perm,
                               UKB_afs$nM1_site0.perm,UKB_afs$nM0_site0.perm)


UKB_afs$s_site1_t_obs<- tot_s(UKB_afs$p1_prime_f.obs, UKB_afs$p1_prime_m.obs, 
                              UKB_afs$p1_sexavg.obs)
UKB_afs$s_site1_t_perm<- tot_s(UKB_afs$p1_prime_f.perm, UKB_afs$p1_prime_m.perm, 
                               UKB_afs$p1_sexavg.perm)
UKB_afs$s_site0_t_obs<- tot_s(UKB_afs$p0_prime_f.obs, UKB_afs$p0_prime_m.obs, 
                              UKB_afs$p0_sexavg.obs)
UKB_afs$s_site0_t_perm<- tot_s(UKB_afs$p0_prime_f.perm, UKB_afs$p0_prime_m.perm, 
                               UKB_afs$p0_sexavg.perm)

UKB_afs$s_site1_f_obs = fec_s(UKB_afs$s_site1_v_obs, UKB_afs$s_site1_t_obs, 
                              UKB_afs$p1_sexavg.obs)
UKB_afs$s_site1_f_perm = fec_s(UKB_afs$s_site1_v_perm, UKB_afs$s_site1_t_perm, 
                               UKB_afs$p1_sexavg.perm)
UKB_afs$s_site0_f_obs = fec_s(UKB_afs$s_site0_v_obs, UKB_afs$s_site0_t_obs, 
                              UKB_afs$p0_sexavg.obs)
UKB_afs$s_site0_f_perm = fec_s(UKB_afs$s_site0_v_perm, UKB_afs$s_site0_t_perm, 
                               UKB_afs$p0_sexavg.perm)


#check if analytical ests agree with optimization results

cb_ests <- UKB_afs %>%
  select(SNP1,
         s_site1_v_obs,s_site1_v_perm,
         s_site1_f_obs,s_site1_f_perm,
         s_site1_t_obs,s_site1_t_perm,
         s_site0_v_obs,s_site0_v_perm,
         s_site0_f_obs,s_site0_f_perm,
         s_site0_t_obs,s_site0_t_perm,
         p1_sexavg.obs,p1_sexavg.perm,
         p0_sexavg.obs,p0_sexavg.perm)
names(cb_ests)[1] <- "Site1_SNP"

cb_test<-merge(cb_ests,UKB_all, by="Site1_SNP")

cor.test(cb_test$s_site1_v_obs,cb_test$s_site1_viability)
cor.test(cb_test$s_site1_v_perm,cb_test$s_site1_viability_perm)
cor.test(cb_test$s_site1_f_obs,cb_test$s_site1_fecundity)
cor.test(cb_test$s_site1_f_perm,cb_test$s_site1_fecundity_perm)


################ Sex-specific fecundity estimates at single sites

## function for sex-specific coefficients

sexsp_fec_s <- function(n1,n0,p) {
  a <- n1*(1-p) - n0*p
  b <- (n0+n1)*(1-p)*p
  s <- a/b
  return(s)
}

# calculate coefficients

UKB_afs$sf_F_site1_obs <- sexsp_fec_s(UKB_afs$nF1_f_site1.obs, UKB_afs$nF0_f_site1.obs,
                                      UKB_afs$p1_female.obs)
UKB_afs$sf_M_site1_obs <- sexsp_fec_s(UKB_afs$nM1_f_site1.obs, UKB_afs$nM0_f_site1.obs,
                                      UKB_afs$p1_male.obs)
UKB_afs$sf_F_site0_obs <- sexsp_fec_s(UKB_afs$nF1_f_site0.obs, UKB_afs$nF0_f_site0.obs,
                                      UKB_afs$p0_female.obs)
UKB_afs$sf_M_site0_obs <- sexsp_fec_s(UKB_afs$nM1_f_site0.obs, UKB_afs$nM0_f_site0.obs,
                                      UKB_afs$p0_male.obs)
UKB_afs$sf_F_site1_perm <- sexsp_fec_s(UKB_afs$nF1_f_site1.perm, UKB_afs$nF0_f_site1.perm,
                                       UKB_afs$p1_female.perm)
UKB_afs$sf_M_site1_perm <- sexsp_fec_s(UKB_afs$nM1_f_site1.perm, UKB_afs$nM0_f_site1.perm,
                                       UKB_afs$p1_male.perm)

UKB_afs$sf_F_site0_perm <- sexsp_fec_s(UKB_afs$nF1_f_site0.perm, UKB_afs$nF0_f_site0.perm,
                                       UKB_afs$p0_female.perm)
UKB_afs$sf_M_site0_perm <- sexsp_fec_s(UKB_afs$nM1_f_site0.perm, UKB_afs$nM0_f_site0.perm,
                                       UKB_afs$p0_male.perm)

#tidy
persite_ests <- UKB_afs %>% select(Chrom, Position=SNP1_Position, SNP=SNP1, MAF=SNP1_MAF,
                   sv = s_site1_v_obs, sv_perm=s_site1_v_perm,
                   sf = s_site1_f_obs, sf_perm=s_site1_f_perm,
                   st = s_site1_t_obs, st_perm=s_site1_t_perm,
                   sf_F = sf_F_site1_obs, sf_M = sf_M_site1_obs,
                   sf_F_perm = sf_F_site1_perm, sf_M_perm=sf_M_site1_perm)

#export
write.table(persite_ests, file = paste0(dr, "UKB_persite_estimates_5-20.txt"),
            sep="\t",
            row.names = FALSE,
            append = FALSE)


#evidence for sexually antagonistic selection

#calculate s*
persite_ests$sf_prod <- persite_ests$sf_F*persite_ests$sf_M
persite_ests$sf_prod_perm <- persite_ests$sf_F_perm*persite_ests$sf_M_perm

#among products with negative signs, count number at top 1%
qtl_exp<- quantile(abs(persite_ests$sf_prod_perm[persite_ests$sf_prod_perm < 0]), 0.99)

t_obs <- sum(abs(persite_ests$sf_prod[persite_ests$sf_prod < 0]) > qtl_exp)
t_perm <-  sum(abs(persite_ests$sf_prod_perm[persite_ests$sf_prod_perm < 0]) > qtl_exp)
b_obs <- sum(abs(persite_ests$sf_prod[persite_ests$sf_prod < 0]) <= qtl_exp)
b_perm <- sum(abs(persite_ests$sf_prod_perm[persite_ests$sf_prod_perm < 0]) <= qtl_exp)
t_counts <- cbind(c(t_obs, t_perm), c(b_obs, b_perm))

#test for excess of antagonistic selection
chisq.test(t_counts)

### Simulations (other allele frequencies)
#Load sim data for three s values

#s=0.1
sim_path_s1_p1 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s1_p1.sim")
sim_path_s1_p2 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s1_p2.sim")
sim_path_s1_p3 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s1_p3.sim")
sim_path_s1_p4 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s1_p4.sim")
sim_path_s1_p5 <- paste0(dr, "simulations/sim_d", 
                          0:5, 
                          "_s1_p5.sim")

#s=0.01
sim_path_s01_p1 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s01_p1.sim")
sim_path_s01_p2 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s01_p2.sim")
sim_path_s01_p3 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s01_p3.sim")
sim_path_s01_p4 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s01_p4.sim")
sim_path_s01_p5 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s01_p5.sim")

#s=0.001
sim_path_s001_p1 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s001_p1.sim")
sim_path_s001_p2 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s001_p2.sim")
sim_path_s001_p3 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s001_p3.sim")
sim_path_s001_p4 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s001_p4.sim")
sim_path_s001_p5 <- paste0(dr, "simulations/sim_d", 
                         0:5, 
                         "_s001_p5.sim")

#load
sims_s1_p1 <- lapply(sim_path_s1_p1, fread, sep="\t", 
                      header=TRUE, check.names=FALSE)
sims_s1_p2 <- lapply(sim_path_s1_p2, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s1_p3 <- lapply(sim_path_s1_p3, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s1_p4 <- lapply(sim_path_s1_p4, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s1_p5 <- lapply(sim_path_s1_p5, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)

sims_s01_p1 <- lapply(sim_path_s01_p1, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s01_p2 <- lapply(sim_path_s01_p2, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s01_p3 <- lapply(sim_path_s01_p3, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s01_p4 <- lapply(sim_path_s01_p4, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s01_p5 <- lapply(sim_path_s01_p5, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)

sims_s001_p1 <- lapply(sim_path_s001_p1, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s001_p2 <- lapply(sim_path_s001_p2, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s001_p3 <- lapply(sim_path_s001_p3, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s001_p4 <- lapply(sim_path_s001_p4, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s001_p5 <- lapply(sim_path_s001_p5, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)

#results tables (output all tables, bias, MSE calculations)
sim_s1_p1_result <- getBias_MSE(sims_s1_p1,st = 0.1,
                                 px_hat = 0.1, 
                                 mse_boots = 1000)
sim_s1_p2_result <- getBias_MSE(sims_s1_p2,st = 0.1,
                                 px_hat = 0.2, 
                                 mse_boots = 1000)
sim_s1_p3_result <- getBias_MSE(sims_s1_p3,st = 0.1,
                                 px_hat = 0.3, 
                                 mse_boots = 1000)
sim_s1_p4_result <- getBias_MSE(sims_s1_p4,st = 0.1,
                                 px_hat = 0.4, 
                                 mse_boots = 1000)
sim_s1_p5_result <- getBias_MSE(sims_s1_p5,st = 0.1,
                                 px_hat = 0.5, 
                                 mse_boots = 1000)

sim_s01_p1_result <- getBias_MSE(sims_s01_p1,st = 0.01,
                                px_hat = 0.1, 
                                mse_boots = 1000)
sim_s01_p2_result <- getBias_MSE(sims_s01_p2,st = 0.01,
                                px_hat = 0.2, 
                                mse_boots = 1000)
sim_s01_p3_result <- getBias_MSE(sims_s01_p3,st = 0.01,
                                px_hat = 0.3, 
                                mse_boots = 1000)
sim_s01_p4_result <- getBias_MSE(sims_s01_p4,st = 0.01,
                                px_hat = 0.4, 
                                mse_boots = 1000)
sim_s01_p5_result <- getBias_MSE(sims_s01_p5,st = 0.01,
                                px_hat = 0.5, 
                                mse_boots = 1000)

sim_s001_p1_result <- getBias_MSE(sims_s001_p1,st = 0.001,
                                px_hat = 0.1, 
                                mse_boots = 1000)
sim_s001_p2_result <- getBias_MSE(sims_s001_p2,st = 0.001,
                                px_hat = 0.2, 
                                mse_boots = 1000)
sim_s001_p3_result <- getBias_MSE(sims_s001_p3,st = 0.001,
                                px_hat = 0.3, 
                                mse_boots = 1000)
sim_s001_p4_result <- getBias_MSE(sims_s001_p4,st = 0.001,
                                px_hat = 0.4, 
                                mse_boots = 1000)
sim_s001_p5_result <- getBias_MSE(sims_s001_p5,st = 0.001,
                                px_hat = 0.5, 
                                mse_boots = 1000)

#export results tables

results_p1<- list(sim_s1_p1_result$Results, 
                 sim_s01_p1_result$Results,
                 sim_s001_p1_result$Results)
bias_p1<- list(sim_s1_p1_result$Bias, 
              sim_s01_p1_result$Bias,
              sim_s001_p1_result$Bias)
mse_p1<- list(sim_s1_p1_result$MSE, 
             sim_s01_p1_result$MSE,
             sim_s001_p1_result$MSE)

results_p2<- list(sim_s1_p2_result$Results, 
                  sim_s01_p2_result$Results,
                  sim_s001_p2_result$Results)
bias_p2<- list(sim_s1_p2_result$Bias, 
               sim_s01_p2_result$Bias,
               sim_s001_p2_result$Bias)
mse_p2<- list(sim_s1_p2_result$MSE, 
              sim_s01_p2_result$MSE,
              sim_s001_p2_result$MSE)

results_p3<- list(sim_s1_p3_result$Results, 
                  sim_s01_p3_result$Results,
                  sim_s001_p3_result$Results)
bias_p3<- list(sim_s1_p3_result$Bias, 
               sim_s01_p3_result$Bias,
               sim_s001_p3_result$Bias)
mse_p3<- list(sim_s1_p3_result$MSE, 
              sim_s01_p3_result$MSE,
              sim_s001_p3_result$MSE)

results_p4<- list(sim_s1_p4_result$Results, 
                  sim_s01_p4_result$Results,
                  sim_s001_p4_result$Results)
bias_p4<- list(sim_s1_p4_result$Bias, 
               sim_s01_p4_result$Bias,
               sim_s001_p4_result$Bias)
mse_p4<- list(sim_s1_p4_result$MSE, 
              sim_s01_p4_result$MSE,
              sim_s001_p4_result$MSE)

results_p5<- list(sim_s1_p5_result$Results, 
                  sim_s01_p5_result$Results,
                  sim_s001_p5_result$Results)
bias_p5<- list(sim_s1_p5_result$Bias, 
               sim_s01_p5_result$Bias,
               sim_s001_p5_result$Bias)
mse_p5<- list(sim_s1_p5_result$MSE, 
              sim_s01_p5_result$MSE,
              sim_s001_p5_result$MSE)