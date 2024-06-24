####### 1. Simulate SAS - interpolation/haplotype
#### updated 4/26/24, JMC

#libraries
library(data.table)

##Import functions
source('/scratch/ukb/scripts/ML_functions.R')

#directory
dr <- "/scratch/ukb/data/"    

####################################################################################
## Functions for simulations

##simulate counts 
sim_counts <- function(num_F,num_M,p1,px,p0,ld,s, 
                                use_r=FALSE, 
                                sampling=FALSE,
                                sample_p=FALSE){
  
  d<-runif(1, min = 0, max = 1) #sample d, uniform between 0-1
  
  if (sample_p == TRUE) {
    
    #sample pair of SNPs
    param_samp <- ukb_r2_mafs[sample(nrow(ukb_r2_mafs), 1), ]
    
    #unknown site frequency
    px <- sample(mafs, 1) #sample from MAF dist
    
    #pairwise MAFs
    p1 <- param_samp$MAF_A  
    p0 <- param_samp$MAF_B
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    #ld from pair
    ld <- sqrt(as.numeric(param_samp$R2)) #in terms of correlation r
    
    D_10 <- D_conv(ld, p1, p0)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    
  } else {
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    D_10 <- ifelse(use_r == TRUE, D_conv(ld, p1, p0), ld)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    
  }
  
  #get hap freqs
  hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k=1)
    
  #save
  f0_0 <- hap_terms$f0_0
  f1_0 <- hap_terms$f1_0
  f2_0 <- hap_terms$f2_0
  f3_0 <- hap_terms$f3_0
  
  f0_1 <- hap_terms$f0_1
  f1_1 <- hap_terms$f1_1
  f2_1 <- hap_terms$f2_1
  f3_1 <- hap_terms$f3_1
  
  fem0 <- ((1-s*px)*f0_0) + ((1+s*qx)*f0_1)
  fem1 <- ((1-s*px)*f1_0) + ((1+s*qx)*f1_1)
  fem2 <- ((1-s*px)*f2_0) + ((1+s*qx)*f2_1)
  fem3 <- ((1-s*px)*f3_0) + ((1+s*qx)*f3_1)
  
  male0 <- ((1+s*px)*f0_0) + ((1-s*qx)*f0_1)
  male1 <- ((1+s*px)*f1_0) + ((1-s*qx)*f1_1)
  male2 <- ((1+s*px)*f2_0) + ((1-s*qx)*f2_1)
  male3 <- ((1+s*px)*f3_0) + ((1-s*qx)*f3_1)
  
  while (any(c(fem0, fem1, fem2, fem3, male0, male1, male2, male3) < 0)) {
    
    #sample pair of SNPs
    param_samp <- ukb_r2_mafs[sample(nrow(ukb_r2_mafs), 1), ]
    
    #unknown site frequency
    px <- sample(mafs, 1) #sample from MAF dist
    
    #pairwise MAFs
    p1 <- param_samp$MAF_A  
    p0 <- param_samp$MAF_B
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    #ld from pair
    ld <- sqrt(as.numeric(param_samp$R2)) #in terms of correlation r
    
    D_10 <- D_conv(ld, p1, p0)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    #hap terms
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k=1)
    
    f0_0 <- hap_terms$f0_0
    f1_0 <- hap_terms$f1_0
    f2_0 <- hap_terms$f2_0
    f3_0 <- hap_terms$f3_0
    
    f0_1 <- hap_terms$f0_1
    f1_1 <- hap_terms$f1_1
    f2_1 <- hap_terms$f2_1
    f3_1 <- hap_terms$f3_1
    
    fem0 <- ((1-s*px)*f0_0) + ((1+s*qx)*f0_1)
    fem1 <- ((1-s*px)*f1_0) + ((1+s*qx)*f1_1)
    fem2 <- ((1-s*px)*f2_0) + ((1+s*qx)*f2_1)
    fem3 <- ((1-s*px)*f3_0) + ((1+s*qx)*f3_1)
    
    male0 <- ((1+s*px)*f0_0) + ((1-s*qx)*f0_1)
    male1 <- ((1+s*px)*f1_0) + ((1-s*qx)*f1_1)
    male2 <- ((1+s*px)*f2_0) + ((1-s*qx)*f2_1)
    male3 <- ((1+s*px)*f3_0) + ((1-s*qx)*f3_1)
    
  }
  
  if (sampling==TRUE){
    
    prob_f <- c(fem0,fem1,fem2,fem3)
    prob_m <- c(male0,male1,male2,male3)
    
    nF <- rmultinom(1, size = num_F, prob = prob_f)
    nM <- rmultinom(1, size = num_M, prob = prob_m)
    
    f_haps <- c(nF[1], nF[2], nF[3], nF[4])
    m_haps <- c(nM[1], nM[2], nM[3], nM[4])
    
  } else {
    
    nF0 <- fem0*num_F
    nF1 <- fem1*num_F
    nF2 <- fem2*num_F
    nF3 <- fem3*num_F
    
    nM0 <- male0*num_M
    nM1 <- male1*num_M
    nM2 <- male2*num_M
    nM3 <- male3*num_M
    
    f_haps <- c(nF0, nF1, nF2, nF3)
    m_haps <- c(nM0, nM1, nM2, nM3)
  }
  
  
  h_counts <- list("females" = f_haps, 
                   "males" = m_haps,
                   "ps" = c(p1,px,p0),
                   "d" = d,
                   "ld_D10" = D_10,
                   "ld_r" = r_conv(D_10,p1,p0),
                   "freqs" = hap_terms)
  
  return(h_counts)
  
}

#Simulate SAS and estimate selection coefficient 

sim_SAS <- function(nf,nm,s,p1,p0,px,ld,d=0.5,use_r=FALSE,
                     sampling=FALSE,
                     sample_p=FALSE){ 
  
  hcnt <- sim_counts(num_F = nf, num_M = nm, s=s,
                              p1 = p1, px = px, p0 = p0,
                              ld=ld,
                              use_r = use_r, 
                              sampling = sampling, sample_p = sample_p)
  
  px_hat = 0.1317675 #median MAF for all imputed SNPs
  
  #counts
  fem_counts <- hcnt$females
  male_counts <- hcnt$males

  #get freqs
  freqs <- hap_freqs(fem_counts, male_counts, p = px_hat, shrinkage = FALSE)
  freqs0 <- freqs$frq_0
  freqs1 <- freqs$frq_1
  
  #Null ML
  ml_null<-lik.function.n.log(s=0,p=px_hat,fem_counts = fem_counts, male_counts = male_counts,
                     freqs0 = freqs0, freqs1 = freqs1)
  
  
  #optimize (dynamic intervals)
  ml <- optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                                  p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                                  freqs0 = freqs0, freqs1 = freqs1)
  
  if (is.nan(ml$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = fem_counts, male_counts = male_counts,
                                                                  freqs0 = freqs0, freqs1 = freqs1))])
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = fem_counts, male_counts = male_counts,
                                                                    freqs0 = freqs0, freqs1 = freqs1))])
    }
    
    ml <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                   p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                   freqs0 = freqs0, freqs1 = freqs1)
  }
  
  #calculate p-values (LRT)
  LRT_stat<- -2 * (as.numeric(ml$objective) - as.numeric(ml_null))
  p_val <- pchisq(abs(LRT_stat), df=1, lower.tail = FALSE)
  
  #out
  out<-list("maximum" = ml$maximum,
            "objective" = ml$objective,
            "null_objective" = ml_null,
            "pvalue" = p_val)
  
  return(out)
}



####################################################################################

###### Load in data

####load MAF and r2 vectors 
mafs <- scan(file = paste0(dr, "UKB_mafs_imputed_filtered.txt")) #distribution of MAFs
r2_vals <- scan(file = paste0(dr, "UKB_r2_genotyped_filtered.txt")) #distribution of r^2 values

ukb_r2_mafs <- fread(paste0(dr, "UKB_mafs_r2_genotyped.txt"), 
                     sep="\t", header = TRUE, 
                     check.names=FALSE) #pairwise table of MAFs and r^2


####Simulate haplotype counts

##parameters
# nf - number of female haplotypes
# nm - number of male haplotypes

# d: #where is unobserved target SNP (in terns of LD)? d=0.5: target in "center" 

# s:  selection coefficient to simulate

# If you want to manually set allele frequencies instead of sampling:
# set "sample_p=FALSE"

# p1: zygote allele frequency at site 1 (leftmost site)
# px: zygote allele frequency of unknown target x
# p0: zygote allele frequency of site 0 (rightmost site)

# ld: input ld (in terms of D) between site 1 and 0. To use r instead of D, set use_r=TRUE. 
#     only use when setting allele frequencies manually

# use_r: TRUE or FALSE, provide ld as correlation between sites 1 and 0 (r) instead of D
#       only set to TRUE when setting allele frequencies manually

# sample_p = TRUE will randomly sample from the data to get ld and the three allele frequencies

# sampling = TRUE will sample from multinomial to get haplotype counts


#TEST (outputs ML and estimate)
simulate1 <- sim_SAS(nf=10000, nm=10000, s=0.01, sample_p = TRUE, sampling=TRUE)

#replicate
sims_100_s <- lapply(1:100, function(x) sim_SAS(nf=10000, nm=10000, s=0.01, sample_p = TRUE)$maximum)
mean(unlist(sims_100_s)) # close to 0.01