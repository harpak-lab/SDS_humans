############# Cole et al, 2024 Plotting
############# updated 4/25/24
#######################################

## packages
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(qqman)
library(scales)
library(ggExtra)

dr <- "/scratch/ukb/data/"

## sim data
sim_path_s1 <- paste0(dr, "simulations/sim_s1_p", 
                         c(1,13,2:5), 
                         ".txt")
sim_path_s01 <- paste0(dr, "simulations/sim_s01_p", 
                      c(1,13,2:5), 
                      ".txt")
sim_path_s001 <- paste0(dr, "simulations/sim_s001_p", 
                      c(1,13,2:5), 
                      ".txt")
sims_s1 <- lapply(sim_path_s1, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s01 <- lapply(sim_path_s01, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)
sims_s001 <- lapply(sim_path_s001, fread, sep="\t", 
                     header=TRUE, check.names=FALSE)

## UKB data
#results
UKB_all <- fread(paste0(dr, "UKB_all_5-22.txt"))
UKB_filt <- fread(paste0(dr, "UKB_filtered_5-22.txt"))
UKB_filt_pruned <- fread(paste0(dr, "UKB_filtered_pruned_5-22.txt")) 
persite_ests <- fread(paste0(dr, "UKB_persite_estimates_5-20.txt"))

## Regression results
GWAS_results <- read.table(paste0(dr, "regression.lfsr.mash.p1.txt"), 
                           header = TRUE)
arm_R_sample<- read.table(paste0(dr, "arm_R_sample.txt", 
                          header=TRUE))
arm_L_sample<- read.table(paste0(dr, "arm_L_sample.txt", 
                          header=TRUE))

## ABC data
abc.df<- read.table(paste0(dr, "ABC_Results_Top1Perc.txt", 
                           header=TRUE))
thresh.ests <-read.table(paste0(dr, "Threshold_vs_ABCestimates.txt", 
                           header=TRUE))
downsampling.ests <- read.table(paste0(dr, "DownSampling_vs_ABCestimates.txt", 
                                                    header=TRUE))

######### FIG 1C, Simulation plot

sim_results_s1 <- sims_s1[[2]]$Results
sim_results_s1$Rel_bias_int <- (sim_results_s1$int_s/0.1)-1
sim_results_s1$Rel_bias_sing <-(sim_results_s1$sing_s/0.1)-1

sim_results_s01 <- sims_s01[[2]]$Results
sim_results_s01$Rel_bias_int <- (sim_results_s01$int_s/0.01)-1
sim_results_s01$Rel_bias_sing <-(sim_results_s01$sing_s/0.01)-1

sim_results_s001 <- sims_s001[[2]]$Results
sim_results_s001$Rel_bias_int <- (sim_results_s001$int_s/0.001)-1
sim_results_s001$Rel_bias_sing <-(sim_results_s001$sing_s/0.001)-1


sim_plot_s1 <- sim_results_s1 %>%
  group_by(d) %>%
  summarise(
    Mean_Int = mean(Rel_bias_int),
    SE_Int = sd(Rel_bias_int) / sqrt(n()),
    Lower_CI_Int = Mean_Int - qt(0.975, df=n()-1) * SE_Int,
    Upper_CI_Int = Mean_Int + qt(0.975, df=n()-1) * SE_Int,
    Mean_Sing = mean(Rel_bias_sing),
    SE_Sing = sd(Rel_bias_sing) / sqrt(n()),
    Lower_CI_Sing = Mean_Sing - qt(0.975, df=n()-1) * SE_Sing,
    Upper_CI_Sing = Mean_Sing + qt(0.975, df=n()-1) * SE_Sing
  )

sim_plot_s1 <-merge(sim_plot_s1,sims_s1[[2]]$MSE, by="d")

sim_plot_s01 <- sim_results_s01 %>%
  group_by(d) %>%
  summarise(
    Mean_Int = mean(Rel_bias_int),
    SE_Int = sd(Rel_bias_int) / sqrt(n()),
    Lower_CI_Int = Mean_Int - qt(0.975, df=n()-1) * SE_Int,
    Upper_CI_Int = Mean_Int + qt(0.975, df=n()-1) * SE_Int,
    Mean_Sing = mean(Rel_bias_sing),
    SE_Sing = sd(Rel_bias_sing) / sqrt(n()),
    Lower_CI_Sing = Mean_Sing - qt(0.975, df=n()-1) * SE_Sing,
    Upper_CI_Sing = Mean_Sing + qt(0.975, df=n()-1) * SE_Sing
  )

sim_plot_s01<-merge(sim_plot_s01,sims_s01[[2]]$MSE, by="d")

sim_plot_s001 <- sim_results_s001 %>%
  group_by(d) %>%
  summarise(
    Mean_Int = mean(Rel_bias_int),
    SE_Int = sd(Rel_bias_int) / sqrt(n()),
    Lower_CI_Int = Mean_Int - qt(0.975, df=n()-1) * SE_Int,
    Upper_CI_Int = Mean_Int + qt(0.975, df=n()-1) * SE_Int,
    Mean_Sing = mean(Rel_bias_sing),
    SE_Sing = sd(Rel_bias_sing) / sqrt(n()),
    Lower_CI_Sing = Mean_Sing - qt(0.975, df=n()-1) * SE_Sing,
    Upper_CI_Sing = Mean_Sing + qt(0.975, df=n()-1) * SE_Sing
  )

sim_plot_s001<-merge(sim_plot_s001,sims_s001[[2]]$MSE, by="d")


# S=0.01
plot(sim_plot_s01$d, 
     sim_plot_s01$Mean_Int, 
     type = "b", lty=1, col = "#008080", 
     xaxt = "n", yaxt="n", 
     ylim = c(-0.4, 0.05),
     xlim=c(0,0.52),
     ylab = "", xlab = "", main = "")
axis(2,at=c(0.0,-0.2,-0.4),labels=c(0.0, -0.2, -0.4), 
     las=2, cex.axis=1.1)
lines(sim_plot_s01$d, sim_plot_s01$Mean_Sing, 
      type = "b", lty=2,  col = "#777777")
arrows(sim_plot_s01$d, sim_plot_s01$Lower_CI_Int, 
       sim_plot_s01$d, sim_plot_s01$Upper_CI_Int, 
       angle = 90, code = 3, length = 0.05, col = "#008080")
arrows(sim_plot_s01$d, sim_plot_s01$Lower_CI_Sing, 
       sim_plot_s01$d, sim_plot_s01$Upper_CI_Sing, 
       angle = 90, code = 3, length = 0.05, col = "#777777")
abline(h = 0, lty = 1)

plot(sim_plot_s01$d, sim_plot_s01$MSE_int, type = "b", 
     lty=1, lwd=1, col = "#008080", xaxt = "n", 
     yaxt = "n", ylim=c(0.000024,0.000045),
     xlim=c(-0.02,0.5))

lines(sim_plot_s01$d, sim_plot_s01$MSE_sing, 
      col = "#777777", type="b", lty=2)
axis(4,at=c(0.00004,0.00003),
     labels=c(expression(4~x~10^-5),expression(3~x~10^-5)), 
     las=2, cex.axis=1.1)
arrows(sim_plot_s1$d, sim_plot_s01$LowerCI_int, 
       sim_plot_s01$d, sim_plot_s01$UpperCI_int, 
       angle = 90, code = 3, length = 0.05, 
       col = "#008080")
arrows(sim_plot_s1$d, sim_plot_s01$LowerCI_sing, 
       sim_plot_s01$d, sim_plot_s01$UpperCI_sing, 
       angle = 90, code = 3, length = 0.05, col = "#777777")

#S=0.001
plot(sim_plot_s001$d, sim_plot_s001$Mean_Int, 
     type = "b", lty=1, col = "#008080", yaxt="n", 
     ylim = c(-0.4, 0.05),
     xlim=c(0,0.52),
     ylab = "", xlab = "", main = "", cex.axis=1.1)
axis(2,at=c(0.0,-0.2,-0.4),labels=c(0.0, -0.2, -0.4), 
     las=2, cex.axis=1.1)
lines(sim_plot_s001$d, sim_plot_s001$Mean_Sing, type = "b", 
      lty=2, col = "#777777")
arrows(sim_plot_s001$d, sim_plot_s001$Lower_CI_Int, 
       sim_plot_s001$d, sim_plot_s001$Upper_CI_Int, 
       angle = 90, code = 3, length = 0.05, col = "#008080")
arrows(sim_plot_s001$d, sim_plot_s001$Lower_CI_Sing, 
       sim_plot_s001$d, sim_plot_s001$Upper_CI_Sing, angle = 90, 
       code = 3, length = 0.05, col = "#777777")
abline(h = 0, lty = 1)




plot(sim_plot_s001$d, sim_plot_s001$MSE_int, type = "b", lty=1, 
     lwd=1, col = "#008080", yaxt="n", 
     ylim=c(0.0000065, 0.000017),
     xlim=c(-0.05,0.5), cex.axis=1.1)
lines(sim_plot_s001$d, sim_plot_s001$MSE_sing, col = "#777777", 
      type="b", lty=2)
axis(4,at=c(0.000016,0.000012,0.000008),
     labels=c(expression(16~x~10^-6),expression(12~x~10^-6),
              expression(8~x~10^-6)), las=2, cex.axis=1.1)

arrows(sim_plot_s001$d, sim_plot_s001$LowerCI_int, 
       sim_plot_s001$d, sim_plot_s001$UpperCI_int, 
       angle = 90, code = 3, length = 0.05, col = "#008080")
arrows(sim_plot_s001$d, sim_plot_s001$LowerCI_sing, 
       sim_plot_s001$d, sim_plot_s001$UpperCI_sing, 
       angle = 90, code = 3, length = 0.05, col = "#777777")



######### FIG 2, Enrichment plot

## function for chi-squared test
perform_chi_squared <- function(df, bin_row_index) {
  bin_obs <- df$obs_count[bin_row_index]
  bin_perm <- df$perm_count[bin_row_index]
  
  total_obs_excl_bin <- sum(df$obs_count) - bin_obs
  total_perm_excl_bin <- sum(df$perm_count) - bin_perm
  
  # Creating the 2x2 contingency table
  matrix <- matrix(c(bin_obs, total_obs_excl_bin, bin_perm, 
                     total_perm_excl_bin), nrow = 2)
  
  # Perform the chi-squared test
  test_result <- chisq.test(matrix)
  
  # Return the p-value
  return(test_result$p.value)
}



##### Viability

#create separate dfs
s_v_observed <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_v)
s_v_permuted <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_v_perm)

s_v_observed$s_v_obs <- abs(s_v_observed$s_v)
s_v_permuted$s_v_perm <- abs(s_v_permuted$s_v_perm)

# binning
bin_edges <- quantile(s_v_permuted$s_v_perm, 
                      probs = seq(0, 1, length.out = 12), na.rm = TRUE)
s_v_permuted$bin <- cut(s_v_permuted$s_v_perm, breaks = bin_edges, 
                        include.lowest = TRUE, labels = FALSE)

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
s_v_observed$bin <- cut(s_v_observed$s_v_obs, breaks = bin_edges, 
                        include.lowest = TRUE, labels = FALSE)

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
via_results_pruned<- merge(obs_counts, perm_counts, by = "bin", all = TRUE)
via_results_pruned$obs_null <- via_results_pruned$obs_count-via_results_pruned$perm_count
via_results_pruned$obs_null_percent <- 
  (via_results_pruned$obs_null / via_results_pruned$obs_count) * 100

# get chi squared p-values
via_results_pruned$p_value <- sapply(1:nrow(via_results_pruned), 
                                     function(x) perform_chi_squared(via_results_pruned, x))
via_results_pruned$bin_label <- round(as.numeric(via_results_pruned$bin)*100, 
                                      digits=2)


##### Fecundity

# create separate dfs
s_f_observed <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_f)
s_f_permuted <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_f_perm)

s_f_observed$s_f_obs <- abs(s_f_observed$s_f)
s_f_permuted$s_f_perm <- abs(s_f_permuted$s_f_perm)


# binning
bin_edges <- quantile(s_f_permuted$s_f_perm, probs = seq(0, 1, length.out = 12), 
                      na.rm = TRUE)
s_f_permuted$bin <- cut(s_f_permuted$s_f_perm, breaks = bin_edges, 
                        include.lowest = TRUE, labels = FALSE)

# bin medians
bin_medians <- s_f_permuted %>%
  group_by(bin) %>%
  summarise(median_val = median(s_f_perm, na.rm = TRUE)) %>%
  mutate(bin_label = as.character(median_val))  # Create bin_label column

bin_medians <- as.data.frame(bin_medians)

# labels
s_f_permuted <- s_f_permuted %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) # Remove intermediate columns

# bin
s_f_observed$bin <- cut(s_f_observed$s_f_obs, breaks = bin_edges, 
                        include.lowest = TRUE, labels = FALSE)

# bin_label to s_t_observed 
s_f_observed <- s_f_observed %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) # Cleanup


obs_counts <- s_f_observed %>%
  group_by(bin) %>%
  summarise(obs_count = n())

obs_counts <- na.omit(obs_counts)

# counts
perm_counts <- s_f_permuted %>%
  group_by(bin) %>%
  summarise(perm_count = n())

# Merge
fec_results_pruned<- merge(obs_counts, perm_counts, 
                           by = "bin", all = TRUE)
fec_results_pruned$obs_null <- fec_results_pruned$obs_count-fec_results_pruned$perm_count
fec_results_pruned$obs_null_percent <- 
  (fec_results_pruned$obs_null / fec_results_pruned$obs_count) * 100

# get chi squared p-values
fec_results_pruned$p_value <- sapply(1:nrow(fec_results_pruned), function(x) perform_chi_squared(fec_results_pruned, x))
fec_results_pruned$bin_label <- round(as.numeric(fec_results_pruned$bin)*100, digits=2)

##### Total

# create separate dfs
s_t_observed <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_t)
s_t_permuted <- UKB_filt_pruned %>% select(Chrom,Site1_SNP,s_t_perm)

s_t_observed$s_t_obs <- abs(s_t_observed$s_t)
s_t_permuted$s_t_perm <- abs(s_t_permuted$s_t_perm)

# binning
bin_edges <- quantile(s_t_permuted$s_t_perm, probs = seq(0, 1, length.out = 12), na.rm = TRUE)
s_t_permuted$bin <- cut(s_t_permuted$s_t_perm, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# bin medians
bin_medians <- s_t_permuted %>%
  group_by(bin) %>%
  summarise(median_val = median(s_t_perm, na.rm = TRUE)) %>%
  mutate(bin_label = as.character(median_val))  # Create bin_label column

bin_medians <- as.data.frame(bin_medians)

# labels
s_t_permuted <- s_t_permuted %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) # Remove intermediate columns

# bin
s_t_observed$bin <- cut(s_t_observed$s_t_obs, breaks = bin_edges, 
                        include.lowest = TRUE, labels = FALSE)

# bin_label to s_t_observed 
s_t_observed <- s_t_observed %>%
  left_join(bin_medians, by = "bin") %>%
  mutate(bin = as.character(bin_label)) %>%
  select(-median_val, -bin_label) # Cleanup


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
                           by = "bin", all = TRUE)
tot_results_pruned$obs_null <- tot_results_pruned$obs_count-tot_results_pruned$perm_count
tot_results_pruned$obs_null_percent <- 
  (tot_results_pruned$obs_null / tot_results_pruned$obs_count) * 100


# get chi squared p-values
tot_results_pruned$p_value <- sapply(1:nrow(tot_results_pruned), 
                                     function(x) perform_chi_squared(tot_results_pruned, x))
tot_results_pruned$bin_label <- round(as.numeric(tot_results_pruned$bin)*100, 
                                      digits=2)



### Final Plot (FIG 2)

via_plot <- ggplot(via_results_pruned, aes(x = factor(bin, labels = bin_label), 
                                           y = obs_null_percent)) +
  geom_point(size = 8, color = "black") + 
  geom_point(size = 7, color = "darkred") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_classic() +
  ylim(-6,6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), # Black border
#        panel.grid.major = element_line(color = "gray", size = 0.5), # Major grid lines
#        panel.grid.minor = element_line(color = "gray", size = 0.25),
        text=element_text(size=20),
plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) 




fec_plot <- ggplot(fec_results_pruned, aes(x = factor(bin, labels = bin_label), 
                                           y = obs_null_percent)) +
  geom_point(size = 8, color = "black") + 
  geom_point(size = 7, color = "gold") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_classic() +
  ylim(-6,6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), # Black border
#        panel.grid.major = element_line(color = "gray", size = 0.5), # Major grid lines
#        panel.grid.minor = element_line(color = "gray", size = 0.25),
        text=element_text(size=20), # Minor grid lines
plot.margin = margin(0.1, 0.1, 0.1, 0.11, "cm")) 



tot_plot <- ggplot(tot_results_pruned, aes(x = factor(bin, labels = bin_label), 
                                           y = obs_null_percent)) +
  geom_point(size = 8, color = "black") + 
  geom_point(size = 7, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_classic() +
  ylim(-6,6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), # Black border
#        panel.grid.major = element_line(color = "gray", size = 0.5), # Major grid lines
#        panel.grid.minor = element_line(color = "gray", size = 0.25),
        text=element_text(size=20),# Minor grid lines
plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) 


#w:644, h:928, pdf: 6.71 x 9.67
plot_grid(via_plot, fec_plot, tot_plot, align = 'v', ncol = 1)

#ggarrange(via_plot, fec_plot, tot_plot, 
#          ncol=1, nrow=3,
#          labels = c("A.", "B.", "C."))


######### FIG 3, Viability vs fecundity

#sturges rule (for binning)
bin_number <- round(1 + 3.322*log(nrow(UKB_filt_pruned)))

#Selection coefficients

# Divide s_f into x quantiles
s_bin <- UKB_filt_pruned %>%
  mutate(bin = ntile(s_f, bin_number)) %>%
  group_by(bin) %>%
  summarize(x_median = median(s_f), y_median = median(s_v), y_sd = sd(s_v))

summary(s_bin$x_median)
summary(s_bin$y_median)

s_bin$x_median <- s_bin$x_median*100
s_bin$y_median <- s_bin$y_median*100

# Calculate breaks for the x-axis labels
n_breaks <- 5
breaks <- seq(1, nrow(s_bin), length.out = n_breaks)
labels <- c("-0.7","-0.09","0","0.09","0.7")

#Selection coefficient plot for viability vs fecundity
vf_s <- ggplot(s_bin, aes(x = as.factor(bin), y = y_median)) +
  geom_point(size=3, fill="#651FFF", pch=21) + 
  geom_smooth(data=s_bin, aes(x=bin), method='lm', color='darkred', se=FALSE) +
  scale_x_discrete(breaks = s_bin$bin[breaks], labels = labels) +
  ylim(-0.02,0.02) +
  theme_classic() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), plot.margin = margin(1, 1, 1, 0.5, "cm"))+
        #panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  theme(plot.title = element_text(size=22, hjust = 0.5)) 



#Z-scores

# Divide Z_f into x quantiles
Z_bin <- UKB_filt_pruned %>%
  mutate(bin = ntile(Z_f, bin_number)) %>%
  group_by(bin) %>%
  summarize(x_median = median(Z_f), y_median = median(Z_v), y_sd = sd(Z_v))

summary(Z_bin$x_median)


# Calculate breaks for the x-axis labels
n_breaks <- 5
breaks <- seq(1, nrow(Z_bin), length.out = n_breaks)
labels <- c("-2.0","-0.7","0","0.7","2.0")



#Zscore plot for viability vs fecundity
vf_z <- ggplot(Z_bin, aes(x = as.factor(bin), y = y_median)) +
  geom_point(size=3, fill="#FE7F00", pch=21) + 
  geom_smooth(data=Z_bin, aes(x=bin), method='lm', color='darkred', se=FALSE) +
  scale_x_discrete(breaks = Z_bin$bin[breaks], labels = labels) +
  theme_classic() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), plot.margin = margin(1, 1, 1, 0.5, "cm"))+
        #panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  theme(plot.title = element_text(size=22, hjust = 0.5)) 

## combine  (w:1000, h:426; 10.42 x 4.54)
plot_grid(vf_z, vf_s, align = 'v', ncol = 2)


######### FIG 4, Regressions

# Subset data by mode
GWAS_viability <- GWAS_results[GWAS_results$Mode == "Viability", ]
GWAS_fecundity <- GWAS_results[GWAS_results$Mode == "Fecundity", ]
GWAS_total <- GWAS_results[GWAS_results$Mode == "Total", ]

# Define factor levels with blank entry at the top
trait_levels <- c("", "","Albumin", "Arm fat-free mass (L)", "Arm fat-free mass (R)", 
                  "BMI", "Calcium", "Creatinine", "Diastolic blood pressure", 
                  "Eosinophil percentage", "Forced vital capacity", "HbA1c", "Height", 
                  "Hip circumference", "IGF1", "Lymphocyte percentage", "Total protein", 
                  "Pulse rate", "Red blood cell count", "SHBG", "Systolic blood pressure", 
                  "Testosterone", "Urate", "Urea", "Waist circumference", "Waist to hip ratio", 
                  "Weight", "Whole body fat mass", "Waist:hip (BMI adjusted)")

top <- data.frame(Trait = "", Mode = NA, mean_slope = NA, SD_slope = NA, Z=NA, 
                  sample_size_regression = NA, SNPs_after_pval_filtering = 0)
top2 <- data.frame(Trait = "", Mode = NA, mean_slope = NA, SD_slope = NA, Z=NA, 
                   sample_size_regression = NA, SNPs_after_pval_filtering = 0)

# Add the empty row to the top of the dataframe
GWAS_viability <- rbind(top, GWAS_viability)
GWAS_fecundity <- rbind(top, GWAS_fecundity)
GWAS_total <- rbind(top, GWAS_total)
GWAS_viability <- rbind(top2, GWAS_viability)
GWAS_fecundity <- rbind(top2, GWAS_fecundity)
GWAS_total <- rbind(top2, GWAS_total)

# Define color vectors
via_color <- "darkred"
fec_color <- "darkgoldenrod1"
tot_color <- "darkgreen"

#trait names
GWAS_viability$Trait2 <- c("", " ","Albumin", "Arm fat-free mass (L)", "Arm fat-free mass (R)", 
                           "BMI", "Calcium", "Creatinine", "Diastolic blood pressure", 
                           "Eosinophil percentage", "Forced vital capacity", "HbA1c", "Height", 
                           "Hip circumference", "IGF1", "Lymphocyte percentage", "Total protein", 
                           "Pulse rate", "Red blood cell count", "SHBG", "Systolic blood pressure",
                           "Testosterone", "Urate", "Urea", "Waist circumference", "Waist to hip ratio",
                           "Weight", "Whole body fat mass", "Waist:hip (BMI adjusted)")
GWAS_fecundity$Trait2 <- c(""," ","Albumin", "Arm fat-free mass (L)", "Arm fat-free mass (R)", 
                           "BMI", "Calcium", "Creatinine", "Diastolic blood pressure", 
                           "Eosinophil percentage", "Forced vital capacity", "HbA1c", "Height", 
                           "Hip circumference", "IGF1", "Lymphocyte percentage", "Total protein", 
                           "Pulse rate", "Red blood cell count", "SHBG", "Systolic blood pressure", 
                           "Testosterone", "Urate", "Urea", "Waist circumference", "Waist to hip ratio", 
                           "Weight", "Whole body fat mass", "Waist:hip (BMI adjusted)")
GWAS_total$Trait2 <- c(""," ","Albumin", "Arm fat-free mass (L)", "Arm fat-free mass (R)", 
                       "BMI", "Calcium", "Creatinine", "Diastolic blood pressure", 
                       "Eosinophil percentage", "Forced vital capacity", "HbA1c", "Height", 
                       "Hip circumference", "IGF1", "Lymphocyte percentage", "Total protein", 
                       "Pulse rate", "Red blood cell count", "SHBG", "Systolic blood pressure", 
                       "Testosterone", "Urate", "Urea", "Waist circumference", "Waist to hip ratio", 
                       "Weight", "Whole body fat mass", "Waist:hip (BMI adjusted)")

# Create ggplots
gw_p1_via <- ggplot(data = GWAS_viability, aes(y = reorder(Trait2, -GWAS_total$Z), x = -Z)) +
  geom_point(size = 4, fill = via_color, pch = 21) +
  labs(title = '', x = '', y = '') +
  xlim(-1.5, 2.3) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', lwd = 1, alpha = .5) +
  geom_vline(xintercept = 1.96, color = 'red', linetype = 'dashed', lwd = 1, alpha = .5) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

gw_p1_fec <- ggplot(data = GWAS_fecundity, aes(y = reorder(Trait2, -GWAS_total$Z), x = -Z)) +
  geom_point(size = 4, fill = fec_color, pch = 21) +
  labs(title = '', x = '', y = '') +
  xlim(-0.8, 2.3) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', lwd = 1, alpha = .5) +
  geom_vline(xintercept = 1.96, color = 'red', linetype = 'dashed', lwd = 1, alpha = .5) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

gw_p1_tot <- ggplot(data = GWAS_total, aes(y = reorder(Trait2, -Z), x = -Z)) +
  geom_point(size = 4, fill = tot_color, pch = 21) +
  labs(title = '', x = '', y = '') +
  xlim(-0.8, 2.3) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', lwd = 1, alpha = .5) +
  geom_vline(xintercept = 1.96, color = 'red', linetype = 'dashed', lwd = 1, alpha = .5) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Combine plots using cowplot
plot_grid(gw_p1_via, plot_grid(gw_p1_fec, gw_p1_tot), rel_widths = c(1, 1.2))


######### FIG 5, ABC heatmap / genetic load

#transform
abc.df$accepted_s <- 10^abc.df$accepted_s
abc.df$accepted_F <- 10^abc.df$accepted_F

#heatmap/ raster
abc_hmap <- ggplot(data = abc.df, aes(x = accepted_s, y = accepted_F)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  theme_classic() + 
  theme(legend.position = 'left') + 
  scale_fill_distiller(direction = 1, palette = "Reds", name="Density",
                       breaks = c(0, 0.04, 0.08, 0.12, 0.16), 
                       labels = c("0", "0.04", "0.08", "0.12", "0.16")) +
  geom_point(data = abc.df, aes(x= accepted_s, y = accepted_F), alpha = 0) +
  geom_point(aes(x = 10^(-3.337279), y = 10^(-0.722032)), size = 4, 
             shape = 21, fill = "white", color = "black") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                breaks = trans_breaks("log10", function(x) 10^x),
                expand=c(0,0)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                breaks = trans_breaks("log10", function(x) 10^x),
                expand=c(0,0)) +
  annotation_logticks() +
  xlab(expression(paste("Selection coefficient, ", italic(s)))) +
  ylab(expression(paste("Fraction of windows, ", italic(F)))) +
  theme(axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

#display marginal distributions
abc_hmap1 <- ggMarginal(abc_hmap , type = "histogram",  
                        color="black", 
                        fill= "#ECECEC")

# 90% Credible intervals
posterior_s <- sort(abc.df$accepted_s)
posterior_F <- sort(abc.df$accepted_F)

n <- length(posterior_s)
ci <- 0.90
n_int <- floor(ci * n)
wd <- rep(NA, n - n_int)

for (i in 1:(n - n_int)) {
  wd[i] <- posterior_s[i + n_int] - posterior_s[i]
}
minInd <- which.min(wd)
c(posterior_s[minInd], posterior_s[minInd + n_int])

for (i in 1:(n - n_int)) {
  wd[i] <- posterior_F[i + n_int] - posterior_F[i]
}
minInd <- which.min(wd)
c(posterior_F[minInd], posterior_F[minInd + n_int])


# Genetic load
# Parameters (from joint posterior )
s <- 10^(-3.337279)
n <- 100000
Fr<- 10^(-0.722032)

# Range of MAFs
p_values <- seq(1e-4, 1e-2, by = 1e-5)

# Calculate genetic load for each allele frequency
load_values <- sapply(p_values, function(p) {
  1 - (1 - s * p)^(n * Fr)
})
loadvp <- data.frame(p_values, load_values)

# Plot
L_sp <- ggplot(loadvp, aes(x = p_values, y = load_values)) +
  geom_line(size = 1, color = "blue") +
  scale_x_log10(breaks = c(1e-5, 1e-4, 1e-3, 1e-2), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks() +
  geom_hline(yintercept = 0.009, lty="dashed") +
  labs(x = "Allele Frequency (p, log scale)", y = "Genetic Load (L)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) 


#Combine heatmap and load
plot_grid(abc_hmap1, L_sp, nrow = 1, ncol = 2, rel_widths = c(1,0.8), 
          rel_heights = c(1,1))








#######################################

## Supplemental plots

#######################################


### B.1 Sims

## Plotting function

sim_plot <- function(results, bias, mse, n, ylims=NULL, px_spec=NULL, 
                     show_legend=FALSE,
                     legend_x=0.5, legend_y=0.5){
  
  sim_results_plot <- pivot_longer(results, 
                                   cols = c("int_s", "sing_s"), 
                                   names_to = "Method", 
                                   values_to = "s")
  st <- unique(results$true_s)
  
  True <- paste0("s = ", st)
  
  ti<- bquote(hat(italic(p))[x]~'='~ .(px_spec)*', '~italic(s)~'='~.(st))
  
  x <- ggplot(sim_results_plot, aes(as.factor(d), s, fill=Method)) + 
    geom_boxplot(outlier.size = 0.1, show.legend = FALSE) + 
    geom_point(size = -1, aes(fill= Method)) +
    scale_fill_manual(values=c("#69b3a2", "grey"), 
                      labels = c("Haplotype", "Single site"))+
    scale_color_manual(values=c("#69b3a2", "grey"), 
                       labels = c("Haplotype", "Single site")) +
    guides(fill = guide_legend(title = NULL, 
                               size=guide_legend(byrow = TRUE),
                               override.aes = list(shape = 15, size=3, 
                                                   color=c("#69b3a2","grey"),
                                                   fill=NA))) +
    labs(title = ti) + 
    xlab("") +
    ylab(expression(Estimated~italic(hat(s)))) +
    geom_hline(aes(yintercept= st), colour= 'red') +
    theme_bw() +
    theme(legend.position = if(show_legend) c(legend_x, legend_y) else "none",
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.text= element_text(size=8),
          legend.spacing.y = unit(0.0001, 'lines'),
          legend.direction = "vertical",
          legend.box.spacing = unit(0.1, "lines"),  # Reduce space between legend items
          legend.spacing = unit(0.1, "lines"),  # Tighter spacing
          legend.key.size = unit(0.5, "lines"),  # Smaller keys
          axis.title.y = element_text(size = 8),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines')) 
  
  if(!is.null(ylims)){
    x <- x + ylim(ylims)
  }
  
  #BIAS
  guides(colour=guide_legend(override.aes=list(shape=17)))
  
  p_bias <- ggplot(bias, aes(x = d), show.legend = FALSE) +
    geom_line(aes(y = RelErr_int, colour = "RelErr_int", 
                  group = "RelErr_int"), show.legend = FALSE) +
    geom_point(aes(y = RelErr_int, colour = "RelErr_int", 
                   group = "RelErr_int"), show.legend = FALSE) +
    geom_line(aes(y = RelErr_sing, colour = "RelErr_sing",
                  group = "RelErr_sing"), show.legend = FALSE) +
    geom_hline(yintercept = 0.0, linetype = 1) +
    geom_point(aes(y = RelErr_sing, color = "RelErr_sing",
                   group = "RelErr_sing"), show.legend = FALSE) +
    scale_colour_manual(values = c("RelErr_int" = "#69b3a2", 
                                   "RelErr_sing" = "grey"),
                        name = "Method") +
    ylab("Relative Bias")+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 8),
          plot.title=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
  
  MSE_plot<- mse %>%
    pivot_longer(
      cols = c(MSE_int, MSE_sing, UpperCI_int, LowerCI_int, 
               UpperCI_sing, LowerCI_sing),
      names_to = c(".value", "Method"),
      names_pattern = "(.*)_(.*)"
    )
  
  p_mse <- ggplot(MSE_plot, aes(x = d, y = MSE, color = Method, 
                                group = Method)) +
    geom_point(aes(shape = Method), size = 3) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
    labs(
      title = "",
      x = expression(True~target~of~selection*','~italic(d)),
      y = "Mean Squared Error"
    ) +
    scale_color_manual(values = c("sing" = "grey", "int" = "#69b3a2"),
                       name = "Method",  
                       labels = c("Interpolation", "Single site")) +
    scale_shape_manual(values = c("sing" = 19, "int" = 19), 
                       name = "Method",  
                       labels = c("Interpolation", "Single site")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 8),
          plot.title=element_blank(),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
  
  plots<-plot_grid(x,p_bias,p_mse,nrow = 3, ncol = 1, align = "v", 
                   rel_heights = c(1,0.8,1))
  plots
  
}

## s=0.1 plots
s1.1 <- sim_plot(results = sims_s1[[1]]$Results, 
                  bias = sims_s1[[1]]$Bias, 
                  mse = sims_s1[[1]]$MSE, n = 300000, 
                  px_spec = 0.1)
s1.13 <- sim_plot(results = sims_s1[[2]]$Results, 
                   bias = sims_s1[[2]]$Bias, 
                  mse = sims_s1[[2]]$MSE, 
                   n = 300000, px_spec = 0.13)
s1.2 <- sim_plot(results = sims_s1[[3]]$Results, 
                  bias = sims_s1[[3]]$Bias, 
                 mse = sims_s1[[3]]$MSE, 
                  n = 300000, px_spec = 0.2)
s1.3 <- sim_plot(results = sims_s1[[4]]$Results, 
                  bias = sims_s1[[4]]$Bias, 
                 mse = sims_s1[[4]]$MSE, 
                  n = 300000, px_spec = 0.3)
s1.4 <- sim_plot(results = sims_s1[[5]]$Results, 
                  bias = sims_s1[[5]]$Bias, 
                 mse = sims_s1[[5]]$MSE, 
                  n = 300000, px_spec = 0.4)
s1.5 <- sim_plot(results = sims_s1[[6]]$Results, 
                  bias = sims_s1[[6]]$Bias, 
                 mse = sims_s1[[6]]$MSE, 
                  n = 300000, px_spec = 0.5)

#w:1060, h:751
ggarrange(s1.1, s1.13, s1.2, s1.3, s1.4, s1.5, ncol=3, nrow=2,
          labels = c("A", "B", "C", "D", "E", "F"))



## s=0.01 plots
s01.1 <- sim_plot(results = sims_s01[[1]]$Results, 
                  bias = sims_s01[[1]]$Bias, 
                  mse = sims_s01[[1]]$MSE, n = 300000, 
                  px_spec = 0.1)
s01.13 <- sim_plot(results = sims_s01[[2]]$Results, 
                   bias = sims_s01[[2]]$Bias, 
                   mse = sims_s01[[2]]$MSE, 
                   n = 300000, px_spec = 0.13)
s01.2 <- sim_plot(results = sims_s01[[3]]$Results, 
                  bias = sims_s01[[3]]$Bias, 
                  mse = sims_s01[[3]]$MSE, 
                  n = 300000, px_spec = 0.2)
s01.3 <- sim_plot(results = sims_s01[[4]]$Results, 
                  bias = sims_s01[[4]]$Bias, 
                  mse = sims_s01[[4]]$MSE, 
                  n = 300000, px_spec = 0.3)
s01.4 <- sim_plot(results = sims_s01[[5]]$Results, 
                  bias = sims_s01[[5]]$Bias, 
                  mse = sims_s01[[5]]$MSE, 
                  n = 300000, px_spec = 0.4)
s01.5 <- sim_plot(results = sims_s01[[6]]$Results, 
                  bias = sims_s01[[6]]$Bias, 
                  mse = sims_s01[[6]]$MSE, 
                  n = 300000, px_spec = 0.5)

#w:1060, h:751
ggarrange(s01.1, s01.13, s01.2, s01.3, s01.4, s01.5, ncol=3, nrow=2,
          labels = c("A", "B", "C", "D", "E", "F"))


## s=0.001 plots
s001.1 <- sim_plot(results = sims_s001[[1]]$Results, 
                  bias = sims_s001[[1]]$Bias, 
                  mse = sims_s001[[1]]$MSE, n = 300000, 
                  px_spec = 0.1)
s001.13 <- sim_plot(results = sims_s001[[2]]$Results, 
                   bias = sims_s001[[2]]$Bias, 
                   mse = sims_s001[[2]]$MSE, 
                   n = 300000, px_spec = 0.13)
s001.2 <- sim_plot(results = sims_s001[[3]]$Results, 
                  bias = sims_s001[[3]]$Bias, 
                  mse = sims_s001[[3]]$MSE, 
                  n = 300000, px_spec = 0.2)
s001.3 <- sim_plot(results = sims_s001[[4]]$Results, 
                  bias = sims_s001[[4]]$Bias, 
                  mse = sims_s001[[4]]$MSE, 
                  n = 300000, px_spec = 0.3)
s001.4 <- sim_plot(results = sims_s001[[5]]$Results, 
                  bias = sims_s001[[5]]$Bias, 
                  mse = sims_s001[[5]]$MSE, 
                  n = 300000, px_spec = 0.4)
s001.5 <- sim_plot(results = sims_s001[[6]]$Results, 
                  bias = sims_s001[[6]]$Bias, 
                  mse = sims_s001[[6]]$MSE, 
                  n = 300000, px_spec = 0.5)


#w:1060, h:751
ggarrange(s001.1, s001.13, s001.2, s001.3, s001.4, s001.5, ncol=3, nrow=2,
          labels = c("A", "B", "C", "D", "E", "F"))


#########

### B.2 Distributions of s hat, empirical
s_v_hist <- ggplot(UKB_filt, aes(x=s_v)) + 
  geom_histogram(fill="darkred", col="black", bins=40) + 
  ggtitle("Viability (observed)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

s_f_hist <- ggplot(UKB_filt, aes(x=s_f)) + 
  geom_histogram(fill="gold", col="black", bins=40) + 
  ggtitle("Fecundity (observed)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15), 
        axis.title.y=element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


s_t_hist <- ggplot(UKB_filt, aes(x=s_t)) + 
  geom_histogram(fill="darkgreen", col="black", bins=40) + 
  ggtitle("Total (observed)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 1, 0.1, 0.1, "cm"))


#(w:1059, h:318)
plot_grid(s_v_hist, s_f_hist, s_t_hist,
          align = 'h', ncol = 3)

### Distributions of s hat, permuted 
s_v_hist_p <- ggplot(UKB_filt, aes(x=s_v_perm)) + 
  geom_histogram(fill="darkred", col="black", bins=40) + 
  ggtitle("Viability (permuted)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

s_f_hist_p <- ggplot(UKB_filt, aes(x=s_f_perm)) + 
  geom_histogram(fill="gold", col="black", bins=40) + 
  ggtitle("Fecundity (permuted)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15), 
        axis.title.y=element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


s_t_hist_p <- ggplot(UKB_filt, aes(x=s_t_perm)) + 
  geom_histogram(fill="darkgreen", col="black", bins=40) + 
  ggtitle("Total (permuted)") +
  xlab("Selection coefficient") +
  ylab("Frequency") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        plot.margin = margin(0.1, 1, 0.1, 0.1, "cm"))


#(w:1059, h:318)
plot_grid(s_v_hist_p, s_f_hist_p, s_t_hist_p,
          align = 'h', ncol = 3)


### B.3 P-value histograms, Manhattan and QQ plots

## histograms, viability

breaks_v <- pretty(range(UKB_filt$s_v_pval), 
              n = nclass.Sturges(UKB_filt$s_v_pval),
              min.n = 1)
via_p_obs <- ggplot(UKB_filt, aes(x=s_v_pval)) + 
  geom_histogram(fill="darkred", col="black", breaks=breaks_v)+
  ggtitle("Viability (observed)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

via_p_perm <- ggplot(UKB_filt, aes(x=s_v_pval_perm)) + 
  geom_histogram(fill="darkred", col="black", breaks=breaks_v)+
  ggtitle("Viability (permuted)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


### QQ plots, viability

#lambda
#lam_v <- median(qchisq(1 - UKB_filt$s_v_pval, df = 1))/median(qchisq(1 - UKB_filt$s_v_pval_perm, df = 1))

#plot
qq_via<-ggplot(data.frame(exp = -log10(sort(UKB_filt$s_v_pval_perm)), 
                          obs = -log10(sort(UKB_filt$s_v_pval)),),
       aes(x = exp, y = obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(x = expression(Permuted~"-"*log[10](italic(p)~value)), 
       y = expression(Observed~"-"*log[10](italic(p)~value)),
       title = "QQ Plot of Observed vs. Permuted (Viability)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
 
 
 
## Combined plot (w:1152, h:1048)
ggarrange(via_p_obs, via_p_perm, qq_via,
          align = 'h', ncol = 2, nrow=2,
          labels = c("A","B","C"))
 

## histograms, fecundity

breaks_f <- pretty(range(UKB_filt$s_f_pval), 
                   n = nclass.Sturges(UKB_filt$s_f_pval),
                   min.n = 1)
fec_p_obs <- ggplot(UKB_filt, aes(x=s_f_pval)) + 
  geom_histogram(fill="gold", col="black", breaks=breaks_f)+
  ggtitle("Fecundity (observed)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

fec_p_perm <- ggplot(UKB_filt, aes(x=s_f_pval_perm)) + 
  geom_histogram(fill="gold", col="black", breaks=breaks_f)+
  ggtitle("Fecundity (permuted)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


### QQ plots, fecundity

#lambda
#lam_f <- median(qchisq(1 - UKB_filt$s_f_pval, df = 1))/median(qchisq(1 - UKB_filt$s_f_pval_perm, df = 1))

#plot
qq_fec<-ggplot(data.frame(exp = -log10(sort(UKB_filt$s_f_pval_perm)), 
                          obs = -log10(sort(UKB_filt$s_f_pval))),
               aes(x = exp, y = obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(x = expression(Permuted~"-"*log[10](italic(p)~value)), 
       y = expression(Observed~"-"*log[10](italic(p)~value)),
       title = "QQ Plot of Observed vs. Permuted (Fecundity)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))



## Combined plot (w:1152, h:1048)
ggarrange(fec_p_obs, fec_p_perm, qq_fec,
          align = 'h', ncol = 2, nrow=2,
          labels = c("A","B","C"))

 
## histogram, total

breaks_t <- pretty(range(UKB_filt$s_t_pval), 
                   n = nclass.Sturges(UKB_filt$s_t_pval),
                   min.n = 1)
tot_p_obs <- ggplot(UKB_filt, aes(x=s_t_pval)) + 
  geom_histogram(fill="darkgreen", col="black", breaks=breaks_t)+
  ggtitle("Total (observed)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

tot_p_perm <- ggplot(UKB_filt, aes(x=s_t_pval_perm)) + 
  geom_histogram(fill="darkgreen", col="black", breaks=breaks_t)+
  ggtitle("Total (permuted)") +
  xlab(expression(italic(p)~values)) +
  ylab("Frequency") + 
  ylim(0,15000) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))


### QQ plots, total

#lambda
#lam_t <- median(qchisq(1 - UKB_filt$s_t_pval, df = 1))/median(qchisq(1 - UKB_filt$s_t_pval_perm, df = 1))

#plot
qq_tot<-ggplot(data.frame(exp = -log10(sort(UKB_filt$s_t_pval_perm)), 
                          obs = -log10(sort(UKB_filt$s_t_pval))),
               aes(x = exp, y = obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(x = expression(Permuted~"-"*log[10](italic(p)~value)), 
       y = expression(Observed~"-"*log[10](italic(p)~value)),
       title = "QQ Plot of Observed vs. Permuted (Total)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=16, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))



## Combined plot (w:1152, h:1048)
ggarrange(tot_p_obs, tot_p_perm, qq_tot,
          align = 'h', ncol = 2, nrow=2,
          labels = c("A","B","C"))


#Manhattan plots


##observed
man_data_v <- UKB_via_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_viability", "SNP" = "Site1_SNP")
man_data_f <- UKB_fec_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_fecundity", "SNP" = "Site1_SNP")
man_data_t <- UKB_tot_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_total", "SNP" = "Site1_SNP")

##bonferonni, 0.05/length(man_data_v$P)

par(mfrow = c(3, 1), mar=c(4.5,4.5,4.5,4.5)) # 
manhattan(man_data_v, cex.axis=1, cex.lab=1.5, cex.main=2, main="Viability", ylim=c(0,16), suggestiveline = FALSE,  
          genomewideline = -log10(0.000000201565))
manhattan(man_data_f, cex.axis=1, cex.lab=1.5, cex.main=2, main="Fecundity", ylim=c(0,16), suggestiveline = FALSE,  
          genomewideline = -log10(0.000000201565))
manhattan(man_data_t, cex.axis=1, cex.lab=1.5, cex.main=2, main="Total", ylim=c(0,16), suggestiveline = FALSE,  
          genomewideline = -log10(0.000000201565))


##permuted
man_data_v_perm <- UKB_via_perm_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_viability", "SNP" = "Site1_SNP")
man_data_f_perm <- UKB_fec_perm_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_fecundity", "SNP" = "Site1_SNP")
man_data_t_perm <- UKB_tot_perm_filt2 %>%
  select("CHR" = "Chrom", "BP" = "Site1_Position", "P"="Pval_int_boot_total", "SNP" = "Site1_SNP")

par(mfrow = c(3, 1), mar=c(4.5,4.5,4.5,4.5)) # 
manhattan(man_data_v_perm, cex.axis=1, cex.lab=1.5, cex.main=2, main="Viability (permuted)", ylim=c(0,16), 
          suggestiveline = FALSE, genomewideline = -log10(0.000000201565))
manhattan(man_data_f_perm , cex.axis=1, cex.lab=1.5, cex.main=2, main="Fecundity (permuted)", ylim=c(0,16), 
          suggestiveline = FALSE, genomewideline = -log10(0.000000201565))
manhattan(man_data_t_perm , cex.axis=1, cex.lab=1.5, cex.main=2, main="Total (permuted)", ylim=c(0,16), 
          suggestiveline = FALSE, genomewideline = -log10(0.000000201565))


## B.5 Evidence for sexually antagonistic selection

#calculate s*
persite_ests$sf_prod <- persite_ests$sf_F*persite_ests$sf_M
persite_ests$sf_prod_perm <- persite_ests$sf_F_perm*persite_ests$sf_M_perm

#bin antagonistic loci (bottom 50)
sf_obs <- persite_ests %>% 
  select(Chrom,SNP,sf_prod) %>%
  filter(sf_prod < 0)
sf_perm <- persite_ests %>% 
  select(Chrom,SNP,sf_prod_perm) %>%
  filter(sf_prod_perm < 0)

sf_obs$sf_prod <- abs(sf_obs$sf_prod)
sf_perm$sf_prod_perm <- abs(sf_perm$sf_prod_perm)

# binning
bin_ed <- quantile(sf_perm$sf_prod_perm, 
                      probs = seq(0, 1, length.out = 101), 
                      na.rm = TRUE)
bin_ed[1] <- 0
bin_ed[101] <- 1

sf_perm$bin <- cut(sf_perm$sf_prod_perm, 
                        breaks = bin_ed, 
                        include.lowest = TRUE, 
                        labels = FALSE)
sf_obs$bin <- cut(sf_obs$sf_prod, 
                        breaks = bin_ed, 
                        include.lowest = TRUE, 
                        labels = FALSE)

sf_obs_counts <- sf_obs %>%
  group_by(bin) %>%
  summarise(obs_count = n())
sf_perm_counts <- sf_perm %>%
  group_by(bin) %>%
  summarise(perm_count = n())
sf_prods <- merge(sf_obs_counts, sf_perm_counts, 
                           by = "bin", 
                           all = TRUE)

sf_prods$obs_null <- sf_prods$obs_count-sf_prods$perm_count
sf_prods$obs_null_percent <- (sf_prods$obs_null / sf_prods$obs_count) * 100

sf_prod_plot<-ggplot(sf_prods, aes(x = bin, y = obs_null_percent)) +
  geom_point(size= 2, color = "black") + # Larger black points for background/border
  geom_point(size = 1, color = "blue") + # Smaller red points on top
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Horizontal line at zero
  theme_classic() +
  geom_smooth(data=subset(sf_prods), method='loess',formula=y~x,se=F,size=1, col="darkred")+
  ggtitle(expression(Quantiles~of~negative~italic(s)[f]^"*"~values)) +
  ylab("Percentage of sites: Observed - Null")+
  xlab(expression(Negative~italic(s)[f]^"*"~(absolute~values)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), # Black border
        text=element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) 

sf_prod_plot
 
## B.5 Compare estimates to Ruzicka et al, 2022
 
#FST DATA (from Ruzicka F, Holman L, Connallon T. Polygenic signals of sex differences in selection in humans from 
#the UK Biobank. PLoS Biol. 2022 Sep 6;20(9):e3001768. 
#doi: 10.1371/journal.pbio.3001768. PMID: 36067235; PMCID: PMC9481184.https://zenodo.org/record/6824671 )

r_fst <- read.table("/stor/work/Kirkpatrick/scratch/Jared/ukb/analysis_new_hap/hap_data/Fst_Data.txt", 
                    sep="", header = TRUE, 
                    check.names=FALSE, 
                    stringsAsFactors = FALSE)
names(r_fst)[1] <- "SNP"

## extract s values for single sites
site1_snps <- UKB_all %>% select(SNP = Site1_SNP, 
                   sv = s_site1_viability, 
                   sf = s_site1_fecundity, 
                   st = s_site1_total)

site0_snps <- UKB_all %>% select(SNP = Site0_SNP, 
                   sv = s_site0_viability, 
                   sf = s_site0_fecundity, 
                   st = s_site0_total)

site_snps<-rbind(site1_snps, site0_snps)
site_snps <- site_snps[!duplicated(site_snps[, 1]), ]

#merge
fst_m <- merge(r_fst, site_snps, by="SNP")

#calculate s coefficients from FST
fst_m$p_hat <- (fst_m$p_M_VIABILITY + fst_m$p_F_VIABILITY)/2
fst_m$p_hat_f <-(fst_m$p_M + fst_m$p_F)/2

fst_m$sr_v <- sqrt(fst_m$FST_VIABILITY/(fst_m$p_hat *(1-fst_m$p_hat)))
fst_m$sr_f <- sqrt(fst_m$FST_REPRODUCTIVE/(fst_m$p_hat_f *(1-fst_m$p_hat_f)))
fst_m$sr_t <- sqrt(fst_m$FST_GAMETIC/(fst_m$p_hat *(1-fst_m$p_hat)))

fst_m$abs_sv <- abs(fst_m$sv)
fst_m$abs_sf <- abs(fst_m$sf)
fst_m$abs_st <- abs(fst_m$st)

#get Pearson's r
cor.test(fst_m$abs_sv,fst_m$sr_v)
cor.test(fst_m$abs_sf,fst_m$sr_f)
cor.test(fst_m$abs_st,fst_m$sr_t)

#Plot
par(mfrow = c(1, 3),
    mar = c(2.5,2.5,2.5,2.5))

plot(abs(fst_m$sv), fst_m$sr_v,
     col=alpha("darkred",0.2), pch=16,
     ylab="", xlab="", xlim = c(0,0.062),
     ylim=c(0,0.062), cex.axis=1.2)

plot(abs(fst_m$sf), fst_m$sr_f,
     col=alpha("#8B8000",0.2), pch=16,
     ylab="", xlab="", xlim = c(0,0.062),
     ylim=c(0,0.062), cex.axis=1.2)

plot(abs(fst_m$st), fst_m$sr_t,
     col=alpha("darkgreen",0.2), pch=16,
     ylab="", xlab="", xlim = c(0,0.062),
     ylim=c(0,0.062), cex.axis=1.2)


# B.6 Trait regression example results (arm fat free mass R and L)

#plot
arm_L <- ggplot(data=arm_L_sample, aes(x=Effect_Diff, y=s_f))+
  geom_point(color='darkblue', alpha = 0.1, size=arm_L_sample$wt*4000)+
  theme_classic() +
  ylim(-0.02,0.02) +
  #xlim(-0.01,0.01) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  geom_abline(intercept = 1.099853e-05, slope = 0.32, color="black", 
              linetype="solid", size=1) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))

arm_R <- ggplot(data=arm_R_sample, aes(x=Effect_Diff, y=s_f))+
  geom_point(color='darkblue', alpha = 0.1, size=arm_R_sample$wt*4000)+
  theme_classic() +
  ylim(-0.02,0.02) +
  #xlim(-0.01,0.01) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  geom_abline(intercept = 4.015922e-05, slope = 0.4, color="black", 
              linetype="solid", size=1) +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))

#w:867, h:338
plot_grid(arm_R, arm_L, nrow=1, ncol=2, align="hv")


# B.7 ABC 

#dists of posteriors
par(mfrow = c(1, 2))
hist(abc.df$accepted_s, col="#ECECEC",
     main=expression(Posterior~distribution~"for"~italic(s)),
     xlab=expression(Accepted~italic(s)~values~(log[10])))
hist(abc.df$accepted_F, col="#ECECEC",
     main=expression(Posterior~distribution~"for"~italic(F)),
     xlab=expression(Accepted~italic(F)~values~(log[10])))

#how threshold affects params
par(mfrow = c(1, 2))
plot(thresh.ests$Acceptance_Threshold,thresh.ests$S_posterior_mean,
     ylim=c(-5,-1), pch=21,  bg=alpha("gray",0.4),
     main=expression("Posterior mean for"~italic(s)~"vs Acceptance Threshold"),
     ylab=expression("Mean posterior"~italic(s)),
     xlab="Acceptance Threshold")
abline(v=0.01, col="red")

plot(thresh.ests$Acceptance_Threshold,thresh.ests$F_posterior_mean,
     ylim=c(-5,-1), pch=21,  bg=alpha("gray",0.4),
     main=expression("Posterior mean for"~italic(F)~"vs Acceptance Threshold"),
     ylab=expression("Mean posterior"~italic(F)),
     xlab="Acceptance Threshold")
abline(v=0.01, col="red")

#Downsampling
downsampling.ests$NumRuns <- downsampling.ests$DownsamplePercent*50000

par(mfrow = c(1, 2))
plot(downsampling.ests$NumRuns,downsampling.ests$S_posterior_mean, 
     pch=21, xlim=c(0,50000),  bg=alpha("gray",0.4),
     main=expression("Posterior mean for"~italic(s)~"vs Number of Runs"),
     ylab=expression("Mean posterior"~italic(s)),
     xlab="Number of runs")

plot(downsampling.ests$NumRuns,downsampling.ests$F_posterior_mean, 
     pch=21, xlim=c(0,50000),  bg=alpha("gray",0.4),
     main=expression("Posterior mean for"~italic(F)~"vs Number of Runs"),
     ylab=expression("Mean posterior"~italic(F)),
     xlab="Number of runs")