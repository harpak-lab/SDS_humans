### Functions used in ML pipeline for UK Biobank data
### updated 1/9/24
### Jared Cole

###############################################

#Get window vectors
getWindows <- function(vec, winsize, overlap) {
  starts = seq(1, length(vec), by=winsize-overlap)
  ends   = starts + winsize - 1
  ends[ends > length(vec)] = length(vec)
  windows <- lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
  names(windows) <- starts
  return(windows)
}

#Get haplotypes
gethaps <- function(wnd, input, s) {
  if (win == 1){
    results = input[input$sex == s, wnd]
  }
  else {
    results = do.call(paste0, input[input$sex == s, wnd], quote=TRUE)
  }
  return(results)
}

#Counting function
counting <- function(x) {
  u <- unique(x);
  data.frame(
    haplotype=u,
    count=sapply(u, function(v) { length(which(x==v)) } )
  )
}

#Count haplotypes (binary to base 10), modified to incorporate LRS
haplocounts <- function(input, lrs, sex) {
  convert = lapply(input, strtoi, base="2")
  counts = Map(function(u, v) 
  {
    # get the frequency count from hap element with table
    # convert to data.frame
    dat <- as.data.frame(table(u))
    # sum of ls element by hap element
    dat$lrs <- rowsum(v, u)[,1]
    # rename the columns
    names(dat)[1:2] <- c("haplotype", "count")
    # return the dataset 
    dat
  }, convert, lrs)
  return(counts)
}

# Get frequency matrix of lrs and haplotype counts
haplrs_matrix <- function(input, lrs, sex) {
  convert = lapply(input, strtoi, base="2")
  counts = Map(function(u, v) 
  {
    # get the frequency count from hap element with table
    # convert to data.frame
    dat <- as.data.frame(u)
    # lrs column
    dat$lrs <- v
    # rename the columns
    names(dat)[1:2] <- c("haplotype", "kids")
    # return the dataset 
    dat
  }, convert, lrs)
  
  mat <- Map(function(x){
    setDT(x)
    m <- suppressMessages(dcast(x, kids~haplotype))
    hapcols <- as.character(unique(x$haplotype))
    m
  },counts)
  
  return(mat)
}

#validates x based on y's haplotypes
haplrs_validate <- function(x,y) {
  val_matrix = Map(function(x,y) {
    #get common haplotype columns
    joined <- suppressMessages(data.table(left_join(x,y)))
    
    #set NAs to 0
    joined[is.na(joined)] <- 0
    
    #order the haplotype columns
    suppressWarnings(setcolorder(joined, c("kids",sort(as.numeric(colnames(joined))))))
    
    #return
    joined
  },x,y)
  
  return(val_matrix)
} 

#Get LRS
getlrs <- function(wnd, input, s) {
  results = as.vector(input[input$sex == s, 3])
  results = as.numeric(unlist(results))
  return(results)
}

#Get snp IDs
getsnppos <- function(wnd, leading, input) {
  as.numeric(unlist(colnames(input[wins[[wnd]][leading]])))
}

# Determine variants in haplotypes
allele <- function(i,j) {
  x = i/(2^(j-1))
  if (floor(x) %% 2 != 0) {
    return(1)
  } else {
    return(0)
  }
}

#Sex averaged frequency, viability
sex_avg_frq <- function(data,window,hap){ 
  win_hap <- data[data$Window == window, ]
  sumF <- sum(win_hap$female_counts)
  sumM <- sum(win_hap$male_counts)
  hapF <- win_hap[win_hap$haplotype == hap, ]$female_counts
  hapM <- win_hap[win_hap$haplotype == hap, ]$male_counts
  hapF <- ifelse(length(hapF) == 0, 0, hapF)
  hapM <- ifelse(length(hapM) == 0, 0, hapM)
  a <- hapF/sumF
  b <- hapM/sumM   
  frq <- (a+b)/2
  return(frq) 
}


#Sex averaged frequency, total
sex_avg_frq2 <- function(data,window,hap){ 
  win_hap <- data[data$Window == window, ]
  sumF <- sum(win_hap$adj_female_counts)
  sumM <- sum(win_hap$adj_male_counts)
  hapF <- win_hap[win_hap$haplotype == hap, ]$adj_female_counts
  hapM <- win_hap[win_hap$haplotype == hap, ]$adj_male_counts
  hapF <- ifelse(length(hapF) == 0, 0, hapF)
  hapM <- ifelse(length(hapM) == 0, 0, hapM)
  a <- hapF/sumF
  b <- hapM/sumM   
  frq <- (a+b)/2
  return(frq) 
}

#function to bootstrap matrix data
haplrs_boot <- function(input) {
  mat <- Map(function(x){
    setDT(x)
    kids <- x$kids
    m<-prop.table(x[,2:ncol(x)])
    m1<-matrix(rmultinom(1,sum(x[,2:ncol(x)]),as.matrix(m)), nrow=nrow(m), ncol=ncol(m))
    rownames(m1) <- rownames(m)
    colnames(m1) <- colnames(m)
    m1 <- cbind(kids,m1)
    m1
  },input)
  return(mat)
}


# D conversion (to get r)
D_conv <- function(r,p1,p2){
  D <- r*(sqrt(p1*(1-p1))*sqrt(p2*(1-p2)))
  return(D)
}

#r conversion (to get D)
r_conv <- function(D,p1,p2){
  r <- D/(sqrt(p1*(1-p1))*sqrt(p2*(1-p2)))
  return(r)
}



# zygote haplotype frequiencies for 3-sites
hap_f_calcs <- function(p0,p1,px,D_1x,D_x0,D_10,k){
  
  q0 <- 1-p0
  q1 <- 1-p1
  qx <- 1-px
  
  #with given value of k
  f0_0 <- q0*q1*qx + qx*D_10 + k*q0*D_1x + k*q1*D_x0
  f1_0 <- p0*q1*qx - qx*D_10 + k*p0*D_1x - k*q1*D_x0
  f2_0 <- p1*q0*qx - qx*D_10 - k*q0*D_1x + k*p1*D_x0
  f3_0 <- p0*p1*qx + qx*D_10 - k*p0*D_1x - k*p1*D_x0
  f0_1 <- px*q0*q1 + px*D_10 - k*q0*D_1x - k*q1*D_x0
  f1_1 <- p0*px*q1 - px*D_10 - k*p0*D_1x + k*q1*D_x0
  f2_1 <- p1*px*q0 - px*D_10 + k*q0*D_1x - k*p1*D_x0
  f3_1 <- p0*p1*px + px*D_10 + k*p0*D_1x + k*p1*D_x0
  
  fq <- list(f0_0 = f0_0, f1_0 = f1_0, f2_0 = f2_0, f3_0 = f3_0, 
             f0_1 = f0_1, f1_1 = f1_1, f2_1 = f2_1, f3_1 = f3_1)
  
  return(fq)
}





### Find allele frequency px (by averageing between two flanking sites, observing 2 sites)

#allele freqs function
pfreq <- function(h){
  
  pF0 <- (h$females[2]+h$females[4])/sum(h$females)
  pF1 <- (h$females[3]+h$females[4])/sum(h$females)
  pFx <- (pF0 + pF1)/2
  
  pM0 <- (h$males[2]+h$males[4])/sum(h$males)
  pM1 <- (h$males[3]+h$males[4])/sum(h$males)
  pMx <- (pM0 + pM1)/2
  
  px <- (pFx + pMx)/2
  
  return(px)
  
}






### Infer the haplotype frequencies in zygotes with 2 observable sites

# get zygote freqs (assumes target in center)
hap_freqs <- function(f,m,p, k=1, shrinkage=FALSE){
  
  #counts
  fcounts <- f
  mcounts <- m
  
  #sums
  NF <- sum(fcounts)
  NM <- sum(mcounts)
  
  #allele freqs
  pF0 <- (fcounts[2]+fcounts[4])/NF
  pF1 <- (fcounts[3]+fcounts[4])/NF
  
  pM0 <- (mcounts[2]+mcounts[4])/NM
  pM1 <- (mcounts[3]+mcounts[4])/NM
  
  qF0 <- 1-pF0
  qF1 <- 1-pF1
  
  qM0 <- 1-pM0
  qM1 <- 1-pM1
  
  p0 <- (pF0+pM0)/2
  p1 <- (pF1+pM1)/2
  px <- p
  
  q0 <- 1-p0
  q1 <- 1-p1
  
  qx <- 1-px
  
  ## LD (infer from two sites, averaged between sexes)
  f4 <- ((fcounts[4]/NF)+(mcounts[4]/NM))/2
  D_10 <- f4 - (p0*p1)
  
  I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
  
  D_1x <- (abs(D_10)^(0.5))*I_D_10*(px*qx)^(0.5)*((p1*q1)/(p0*q0))^(1/4)
  D_x0 <- (abs(D_10)^(0.5))*(px*qx)^(0.5)*((p0*q0)/(p1*q1))^(1/4)
  
  r10 <- D_10/sqrt(p0*q0*p1*q1)
  r1x <- D_1x/sqrt(p1*q1*px*qx)
  rx0 <- D_x0/sqrt(p0*q0*px*qx)
  
  #starting haplotype freqs, determine k
  
  if (shrinkage==TRUE){
    
    k<-1
    
    repeat {
      hap_terms <- hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
      
      k_values <- numeric() # store k
      
      for (term_name in names(hap_terms)) {
        term_value <- hap_terms[[term_name]]
        if (!is.na(term_value) && term_value < 0) {
          if (term_name == "f0_0") {
            k_values <- c(k_values, (-q0*q1*qx - qx*D_10)/(q0*D_1x + q1*D_x0))
          } else if (term_name == "f1_0") {
            k_values <- c(k_values, (-p0*q1*qx + qx*D_10)/(p0*D_1x - q1*D_x0))
          } else if (term_name == "f2_0") {
            k_values <- c(k_values, (-p1*q0*qx + qx*D_10)/(-q0*D_1x + p1*D_x0))
          } else if (term_name == "f3_0") {
            k_values <- c(k_values, (p0*p1*qx + qx*D_10)/(p0*D_1x + p1*D_x0))
          } else if (term_name == "f0_1") {
            k_values <- c(k_values, (px*q0*q1 + px*D_10)/(q0*D_1x + q1*D_x0))
          } else if (term_name == "f1_1") {
            k_values <- c(k_values, (p0*px*q1 - px*D_10)/(p0*D_1x - q1*D_x0))
          } else if (term_name == "f2_1") {
            k_values <- c(k_values, (p1*px*q0 - px*D_10)/(-q0*D_1x + p1*D_x0))
          } else if (term_name == "f3_1") {
            k_values <- c(k_values, (-p0*p1*px - px*D_10)/(p0*D_1x + p1*D_x0))
          }
        }
      }
      
      # If all terms are non-negative, we exit the loop
      if (!any(is.na(unlist(hap_terms))) && all(unlist(hap_terms) >= 0)) {
        break
      } else {
        # Update k to the minimum positive k value calculated
        k <- min(k_values[k_values > 0])
        k <- as.numeric(substr(as.character(k), 1, 8))
      }
    }
    
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
    
  } else {
    
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
    
  }
  
  hapfreqs_0 <- c(hap_terms$f0_0,hap_terms$f1_0,hap_terms$f2_0,hap_terms$f3_0)
  hapfreqs_1 <- c(hap_terms$f0_1,hap_terms$f1_1,hap_terms$f2_1,hap_terms$f3_1)
  rvals <- c(r10,r1x,rx0)
  Dvals <- c(D_10,D_1x,D_x0)
  
  frqs <- list("frq_0" = hapfreqs_0,
               "frq_1" = hapfreqs_1,
               "rvals" = rvals,
               "Dvals" = Dvals)
  
  return(frqs)
}



lik.function.n.log <- function(s, p, fem_counts, male_counts, freqs0, freqs1){
  
  #hap freqs
  f0_0 <- freqs0[1]
  f1_0 <- freqs0[2]
  f2_0 <- freqs0[3]
  f3_0 <- freqs0[4]
  f0_1 <- freqs1[1]
  f1_1 <- freqs1[2]
  f2_1 <- freqs1[3]
  f3_1 <- freqs1[4]
  
  #hap counts
  nF0 <- fem_counts[1]
  nF1 <- fem_counts[2]
  nF2 <- fem_counts[3]
  nF3 <- fem_counts[4]
  
  nM0 <- male_counts[1]
  nM1 <- male_counts[2]
  nM2 <- male_counts[3]
  nM3 <- male_counts[4]
  
  #alleles
  px <- p
  qx <- 1-px  
  
  
  #after selection (take logs)
  F0 <- ifelse(nF0 == 0, 0, log(((1-s*px)*f0_0)+((1+s*qx)*f0_1))*nF0)
  F1 <- ifelse(nF1 == 0, 0, log(((1-s*px)*f1_0)+((1+s*qx)*f1_1))*nF1)
  F2 <- ifelse(nF2 == 0, 0, log(((1-s*px)*f2_0)+((1+s*qx)*f2_1))*nF2)
  F3 <- ifelse(nF3 == 0, 0, log(((1-s*px)*f3_0)+((1+s*qx)*f3_1))*nF3)
  
  M0 <- ifelse(nM0 == 0, 0, log(((1+s*px)*f0_0)+((1-s*qx)*f0_1))*nM0)
  M1 <- ifelse(nM1 == 0, 0, log(((1+s*px)*f1_0)+((1-s*qx)*f1_1))*nM1)
  M2 <- ifelse(nM2 == 0, 0, log(((1+s*px)*f2_0)+((1-s*qx)*f2_1))*nM2)
  M3 <- ifelse(nM3 == 0, 0, log(((1+s*px)*f3_0)+((1-s*qx)*f3_1))*nM3)
  
  #log likelihood
  LL <- sum(F0,F1,F2,F3,M0,M1,M2,M3)
  
  return(LL)
}












#bootstrap counts data

## viability 
bootstrap.v <- function(data, pxhat){
  hap_boot <- data %>% group_by(Window) %>% 
    mutate(fM = male_counts/sum(male_counts)) %>%
    mutate(male_counts = as.vector(rmultinom(1, sum(male_counts), fM))) %>%
    mutate(fF = female_counts/sum(female_counts)) %>%
    mutate(female_counts = as.vector(rmultinom(1, sum(female_counts), fF))) %>% 
    select(Chrom,Window,haplotype,female_counts,male_counts,adj_female_counts,adj_male_counts,
           Leading_Position, Leading_SNP, Leading_MAF, Leading_Allele0, Leading_Allele1)
  hap_boot[order(hap_boot$Window, hap_boot$haplotype),] -> hap_boot
  setDT(hap_boot)
  
  winlist <- as.character(unique(hap_boot$Window))
  boot_ml_result <- vector("list", length=length(winlist))
  boot_ml_result <- lapply(winlist, function(window) likelihood.out.v(hap_boot, window, pxhat = pxhat, shrinkage=FALSE))
  
  boot_ml_output <- do.call(rbind, boot_ml_result)
  boot_ml_output[order(boot_ml_output$Window),] -> boot_ml_output
  boot_ml_output <- boot_ml_output %>% 
    select(Window, Site1_Position=Leading_Position, Site1_SNP=Leading_SNP, s_v_int_boot=s_v_int)
  
  return(boot_ml_output)
}

## total
bootstrap.t <- function(data,m,f,pxhat){
  
  #bootstrap matrix
  boot_males <- haplrs_boot(m)
  boot_females <- haplrs_boot(f)
  
  #create counts table
  hap_boot <- data %>% group_by(Window) %>% 
    arrange(Window, haplotype) %>%
    mutate(male_lrs_2 = as.vector(colSums(boot_males[[as.numeric(unique(Window))]][,1] * boot_males[[as.numeric(unique(Window))]][,2:ncol(boot_males[[as.numeric(unique(Window))]])]))) %>%
    mutate(female_lrs_2 = as.vector(colSums(boot_females[[as.numeric(unique(Window))]][,1] * boot_females[[as.numeric(unique(Window))]][,2:ncol(boot_females[[as.numeric(unique(Window))]])])))
  
  #project haplotype counts from lrs
  hap_boot$female_lrs2_2 <- hap_boot$female_lrs_2/2 #halve the total LRS
  hap_boot$male_lrs2_2 <- hap_boot$male_lrs_2/2
  hap_boot <- hap_boot %>% group_by(Window) %>% #get avg LRS for both sexes
    mutate(aF2 = sum(female_lrs2_2)/sum(female_counts))
  hap_boot <- hap_boot %>% group_by(Window) %>% 
    mutate(aM2 = sum(male_lrs2_2)/sum(male_counts))
  hap_boot$adj_female_counts <- hap_boot$female_lrs2_2/(hap_boot$aF2) #weight haplotype
  hap_boot$adj_male_counts <- hap_boot$male_lrs2_2/(hap_boot$aM2)
  
  #save columns
  hap_boot[order(hap_boot$Window, hap_boot$haplotype),] -> hap_boot
  setDT(hap_boot)
  
  winlist <- as.character(unique(hap_boot$Window))
  
  boot_ml_result <- vector("list", length=length(winlist))
  boot_ml_result <- lapply(winlist, function(window) likelihood.out.t(hap_boot, window, 
                                                                      pxhat = pxhat, shrinkage=FALSE))
  boot_ml_output <- do.call(rbind, boot_ml_result)
  boot_ml_output[order(boot_ml_output$Window),] -> boot_ml_output
  boot_ml_output <- boot_ml_output %>% 
    select(Window, Site1_Position=Leading_Position, Site1_SNP=Leading_SNP, s_t_int_boot=s_t_int)
  
  return(boot_ml_output)
}



## fecundity
bootstrap.f <- function(data,m,f,pxhat){
  
  #bootstrap matrix
  boot_males <- haplrs_boot(m)
  boot_females <- haplrs_boot(f)
  
  #create counts table
  hap_boot <- data %>% group_by(Window) %>% 
    arrange(Window, haplotype) %>%
    mutate(male_counts_2 = as.vector(colSums(boot_males[[unique(Window)]][,2:ncol(boot_males[[unique(Window)]])]))) %>%
    mutate(female_counts_2 = as.vector(colSums(boot_females[[as.numeric(unique(Window))]][,2:ncol(boot_females[[as.numeric(unique(Window))]])]))) %>%
    mutate(male_lrs_2 = as.vector(colSums(boot_males[[as.numeric(unique(Window))]][,1] * boot_males[[as.numeric(unique(Window))]][,2:ncol(boot_males[[as.numeric(unique(Window))]])]))) %>%
    mutate(female_lrs_2 = as.vector(colSums(boot_females[[as.numeric(unique(Window))]][,1] * boot_females[[as.numeric(unique(Window))]][,2:ncol(boot_females[[as.numeric(unique(Window))]])])))
  
  #project haplotype counts from lrs
  hap_boot$female_lrs2_2 <- hap_boot$female_lrs_2/2 #halve the total LRS
  hap_boot$male_lrs2_2 <- hap_boot$male_lrs_2/2
  hap_boot <- hap_boot %>% group_by(Window) %>% #get avg LRS for both sexes
    mutate(aF2 = sum(female_lrs2_2)/sum(female_counts))
  hap_boot <- hap_boot %>% group_by(Window) %>% 
    mutate(aM2 = sum(male_lrs2_2)/sum(male_counts))
  hap_boot$adj_female_counts <- hap_boot$female_lrs2_2/(hap_boot$aF2) #weight haplotype
  hap_boot$adj_male_counts <- hap_boot$male_lrs2_2/(hap_boot$aM2)
  
  #save columns
  hap_boot$male_counts <- hap_boot$male_counts_2
  hap_boot$female_counts <- hap_boot$female_counts_2
  hap_boot[order(hap_boot$Window, hap_boot$haplotype),] -> hap_boot
  setDT(hap_boot)
  
  winlist <- as.character(unique(hap_boot$Window))
  
  boot_ml_result <- vector("list", length=length(winlist))
  boot_ml_result <- lapply(winlist, function(window) likelihood.out.f(hap_boot, window, 
                                                                      pxhat = pxhat, shrinkage=FALSE))
  boot_ml_output <- do.call(rbind, boot_ml_result)
  boot_ml_output[order(boot_ml_output$Window),] -> boot_ml_output
  boot_ml_output <- boot_ml_output %>% 
    select(Window, Site1_Position=Leading_Position, Site1_SNP=Leading_SNP, s_f_int_boot=s_f_int)
  
  return(boot_ml_output)
}


#p-values and SEs from bootstrapping

## ses
sderr <- function(est,dist){
  se <- sqrt(sum((dist-est)^2/(length(dist)-1)))
  return(se)
}

## viability

get_boot_p.v <- function(snpid,bootdist,est){
  #point estimates
  s_hat_v <- as.numeric(est[est$Site1_SNP == snpid,]$s_v_int)
  
  #bootstrap distributions
  s_dist_v <- as.vector(as.numeric(bootdist[bootdist$Site1_SNP == snpid, ]$s_v_int_boot))
  
  #SEs
  SE_v <- sd(s_dist_v)
  #SE_v <- sderr(beta_hat_v,betas_dist_v)
  
  #Z-scores
  Z_v <-  s_hat_v/sd(s_dist_v)
  
  #p-values
  P_v <- 2 * pnorm(abs(Z_v),mean = 0, sd = 1, lower.tail = FALSE)
  
  #s_hat falls outside boot dist?
  outside_check <- s_hat_v < min(s_dist_v) | s_hat_v > max(s_dist_v)
  outside <- ifelse(outside_check, "Yes", "No")
  
  out<- data.frame("Site1_SNP" = snpid,
                   "s_int_hat_v" = s_hat_v,
                   "SE_v" = SE_v,
                   "Z_v" = Z_v,
                   "P_v" = P_v,
                   "Outside" = outside)
  return(out)
}

## total
get_boot_p.t <- function(snpid,bootdist,est){
  #point estimates
  s_hat_t <- as.numeric(est[est$Site1_SNP == snpid,]$s_t_int)
  
  #bootstrap distributions
  s_dist_t <- as.vector(as.numeric(bootdist[bootdist$Site1_SNP == snpid, ]$s_t_int_boot))
  
  #SEs
  SE_t <- sd(s_dist_t)
  #SE_t <- sderr(beta_hat_t,betas_dist_t)
  
  #Z-scores
  Z_t <-  s_hat_t/sd(s_dist_t)
  
  #p-values
  P_t <- 2 * pnorm(abs(Z_t),mean = 0, sd = 1, lower.tail = FALSE)
  
  #s_hat falls outside boot dist?
  outside_check <- s_hat_t < min(s_dist_t) | s_hat_t > max(s_dist_t)
  outside <- ifelse(outside_check, "Yes", "No")
  
  out<- data.frame("Site1_SNP" = snpid,
                   "s_int_hat_t" = s_hat_t,
                   "SE_t" = SE_t,
                   "Z_t" = Z_t,
                   "P_t" = P_t,
                   "Outside" = outside)
  return(out)
}


## fecundity
get_boot_p.f <- function(snpid,bootdist,est){
  #point estimates
  s_hat_f <- as.numeric(est[est$Site1_SNP == snpid,]$s_f_int)
  
  #bootstrap distributions
  s_dist_f <- as.vector(as.numeric(bootdist[bootdist$Site1_SNP == snpid, ]$s_f_int_boot))
  
  #SEs
  SE_f <- sd(s_dist_f)
  #SE_f <- sderr(beta_hat_f,betas_dist_f)
  
  #Z-scores
  Z_f <-  s_hat_f/sd(s_dist_f)
  
  #p-values
  P_f <- 2 * pnorm(abs(Z_f),mean = 0, sd = 1, lower.tail = FALSE)
  
  #s_hat falls outside boot dist?
  outside_check <- s_hat_f < min(s_dist_f) | s_hat_f > max(s_dist_f)
  outside <- ifelse(outside_check, "Yes", "No")
  
  out<- data.frame("Site1_SNP" = snpid,
                   "s_int_hat_f" = s_hat_f,
                   "SE_f" = SE_f,
                   "Z_f" = Z_f,
                   "P_f" = P_f,
                   "Outside" = outside)
  return(out)
}


#Negative likelihood function (when minimizing to get maximum)
neg_lik_function <- function(s, p, fem_counts, male_counts, freqs0, freqs1){
  -lik.function.n.log(s, p, fem_counts, male_counts, freqs0, freqs1)
}

# Likelihood function 
lik.function <- function(B, params, freqs){ # n1=nF0, n2= nM0, n3= nF1, n4 = nM1, c = c
  n1 <- params[1]
  n2 <- params[2]
  n3 <- params[3]
  n4 <- params[4]
  c <- params[5]
  p <- freqs[1]     # freqs = p, q
  q <- freqs[2]
  logl <- (n1 * log(1-p*B)) + (n2 * log(1+p*B)) +
    (n3 * log(1+q*B)) + (n4 * log(1-q*B)) + c
  return(logl)
}

# Likelihood function (single)
lik.function.single <- function(s, counts){ # n1=nF0, n2= nM0, n3= nF1, n4 = nM1, c = c
  nF0 <- counts[1]
  nM0 <- counts[2]
  nF1 <- counts[3]
  nM1 <- counts[4]
  pF <- nF1 / (nF0+nF1)
  pM <- nM1 / (nM0+nM1)
  p <- (pM+pF)/2
  q <- 1-p
  logl <- (nF0 * log(1-p*s)) + (nM0 * log(1+p*s)) +
    (nF1 * log(1+q*s)) + (nM1 * log(1-q*s)) + 
    (nF1+nM1)*log(p) + (nF0+nM0)*log(q)
  return(logl)
}

#Optimization function for interpolation (viability)
likelihood.out.v <- function(data,window, pxhat, shrinkage=FALSE){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- as.numeric(unique(unlist(data[data$Window == window,3])))
  expected_haps <- c(0,1,2,3)
  
  #variables 
  fem_counts <- NULL
  male_counts <- NULL
  freqs0 <- NULL
  freqs1 <- NULL
  px_hat <- pxhat
  
  #get counts
  window_counts <- data[data$Window == window,3:5]
  missing_haplotypes <- setdiff(expected_haps, window_counts$haplotype)
  missing_data <- data.frame(haplotype = missing_haplotypes,
                             female_counts = rep(0, length(missing_haplotypes)),
                             male_counts = rep(0, length(missing_haplotypes)))
  combined_data <- rbind(window_counts, missing_data)
  
  fem_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "female_counts"]))
  male_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "male_counts"]))
  
  #get freqs
  freqs <- hap_freqs(fem_counts, male_counts, p = px_hat, shrinkage = shrinkage)
  freqs0 <- freqs$frq_0
  freqs1 <- freqs$frq_1
  
  #optimize (dynamic intervals)
  ml <- suppressWarnings(optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                 p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                 freqs0 = freqs0, freqs1 = freqs1))

  if (is.nan(ml$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = fem_counts, male_counts = male_counts,
                                                                  freqs0 = freqs0, freqs1 = freqs1))]))
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = fem_counts, male_counts = male_counts,
                                                                    freqs0 = freqs0, freqs1 = freqs1))]))
    }

    ml <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                   p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                   freqs0 = freqs0, freqs1 = freqs1)
  }


  #site 1 (leftmost) optimize (single)
  scounts1 <- c(fem_counts[1]+fem_counts[2],
                male_counts[1]+male_counts[2],
                fem_counts[3]+fem_counts[4],
                male_counts[3]+male_counts[4])
  
  ml_st1 <- optimize(lik.function.single,
                     scounts1,
                     interval = c(-0.5,0.5), 
                     maximum = TRUE)
  
  #site 0 (rightmost) optimize (single)
  
  scounts0 <- c(fem_counts[1]+fem_counts[3],
                male_counts[1]+male_counts[3],
                fem_counts[2]+fem_counts[4],
                male_counts[2]+male_counts[4])
  
  ml_st0 <- optimize(lik.function.single,
                     scounts0,
                     interval = c(-0.5,0.5), 
                     maximum = TRUE)
  
  
  ### get results
  s_int <- ml$maximum
  maxl_int <- ml$objective
  s_st1 <- ml_st1$maximum
  maxl_st1 <- ml_st1$objective
  s_st0 <- ml_st0$maximum
  maxl_st0 <- ml_st0$objective
  
  ### null likelihoods
  null_ML_int <- lik.function.n.log(s=0, p = px_hat, 
                                    fem_counts, male_counts, freqs0, freqs1)
  
  null_ML_st1 <- lik.function.single(s=0, scounts1)
  null_ML_st0 <- lik.function.single(s=0, scounts0)
  
  #calculate p-values (Likelihood ratio tests)
  LRT_stat_int <- -2 * (as.numeric(maxl_int) - as.numeric(null_ML_int))
  LRT_stat_st1 <- -2 * (as.numeric(maxl_st1) - as.numeric(null_ML_st1))
  LRT_stat_st0 <- -2 * (as.numeric(maxl_st0) - as.numeric(null_ML_st0))
  
  pval_int <- pchisq(abs(LRT_stat_int), df=1, lower.tail = FALSE)
  pval_st1 <- pchisq(abs(LRT_stat_st1), df=1, lower.tail = FALSE)
  pval_st0 <- pchisq(abs(LRT_stat_st0), df=1, lower.tail = FALSE)
  
  #output info
  chr <- unique(data$Chrom)
  pos <- unique(data[data$Window == window, ]$Leading_Position)
  snpid <- unique(data[data$Window == window, ]$Leading_SNP)
  lmaf <- unique(data[data$Window == window, ]$Leading_MAF)
  al0 <- unique(data[data$Window == window, ]$Leading_Allele0)
  al1 <- unique(data[data$Window == window, ]$Leading_Allele1)
  r10 <- freqs$rvals[1]
  r1x <- freqs$rvals[2]
  rx0 <- freqs$rvals[3]
  D10 <- freqs$Dvals[1]
  D1x <- freqs$Dvals[2]
  Dx0 <- freqs$Dvals[3]
  
  #Output results 
  out <- data.frame(chr,window,pos,snpid,lmaf,al0,al1,r10,r1x,rx0,D10,D1x,Dx0,   
                    s_int, maxl_int, s_st1, maxl_st1, s_st0, maxl_st0,
                    pval_int, pval_st1, pval_st0)
  
  colnames(out) <- (c("Chrom","Window","Leading_Position","Leading_SNP","Leading_MAF","Leading_Allele0","Leading_Allele1",
                      "r_10", "r_1x", "r_x0", "D_10","D_1x","D_x0",
                      "s_v_int","ML_v_int","s_v_site1","ML_v_site1","s_v_site0","ML_v_site0",
                      "Pval_v_int_lrt", "Pval_v_site1_lrt", "Pval_v_site0_lrt"))
  
  output <- rbind(output, out)
  
  return(output)
}

#Optimization function for interpolation (total)
likelihood.out.t<- function(data,window, pxhat, shrinkage=FALSE){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- as.numeric(unique(unlist(data[data$Window == window,3])))
  expected_haps <- c(0,1,2,3)
  
  #variables 
  adj_fem_counts <- NULL
  adj_male_counts <- NULL
  freqs0 <- NULL
  freqs1 <- NULL
  px_hat <- pxhat
  
  #get counts
  window_counts <- data[data$Window == window,c(3,10:11)]
  missing_haplotypes <- setdiff(expected_haps, window_counts$haplotype)
  missing_data <- data.frame(haplotype = missing_haplotypes,
                             adj_female_counts = rep(0, length(missing_haplotypes)),
                             adj_male_counts = rep(0, length(missing_haplotypes)))
  combined_data <- rbind(window_counts, missing_data)
  
  adj_fem_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "adj_female_counts"]))
  adj_male_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "adj_male_counts"]))
  
  #get freqs
  freqs <- hap_freqs(adj_fem_counts, adj_male_counts, p = px_hat, shrinkage = shrinkage)
  freqs0 <- freqs$frq_0
  freqs1 <- freqs$frq_1
  
  #optimize (dynamic intervals)
  ml <- suppressWarnings(optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                                  p = px_hat, fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                  freqs0 = freqs0, freqs1 = freqs1))
  
  if (is.nan(ml$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                                                  freqs0 = freqs0, freqs1 = freqs1))]))
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                                                    freqs0 = freqs0, freqs1 = freqs1))]))
    }
    
    ml <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                   p = px_hat, fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                   freqs0 = freqs0, freqs1 = freqs1)
  }
  
  
  #site 1 (leftmost) optimize (single)
  scounts1 <- c(adj_fem_counts[1]+adj_fem_counts[2],
                adj_male_counts[1]+adj_male_counts[2],
                adj_fem_counts[3]+adj_fem_counts[4],
                adj_male_counts[3]+adj_male_counts[4])
  
  ml_st1 <- optimize(lik.function.single,
                     scounts1,
                     interval = c(-0.5,0.5), 
                     maximum = TRUE)
  
  #site 0 (rightmost) optimize (single)
  
  scounts0 <- c(adj_fem_counts[1]+adj_fem_counts[3],
                adj_male_counts[1]+adj_male_counts[3],
                adj_fem_counts[2]+adj_fem_counts[4],
                adj_male_counts[2]+adj_male_counts[4])
  
  ml_st0 <- optimize(lik.function.single,
                     scounts0,
                     interval = c(-0.5,0.5), 
                     maximum = TRUE)
  
  
  ### get results
  s_int <- ml$maximum
  maxl_int <- ml$objective
  s_st1 <- ml_st1$maximum
  maxl_st1 <- ml_st1$objective
  s_st0 <- ml_st0$maximum
  maxl_st0 <- ml_st0$objective
  
  ### null likelihoods
  null_ML_int <- lik.function.n.log(s=0, p = px_hat, 
                                    adj_fem_counts, adj_male_counts, freqs0, freqs1)
  
  null_ML_st1 <- lik.function.single(s=0, scounts1)
  null_ML_st0 <- lik.function.single(s=0, scounts0)
  
  #calculate p-values (Likelihood ratio tests)
  LRT_stat_int <- -2 * (as.numeric(maxl_int) - as.numeric(null_ML_int))
  LRT_stat_st1 <- -2 * (as.numeric(maxl_st1) - as.numeric(null_ML_st1))
  LRT_stat_st0 <- -2 * (as.numeric(maxl_st0) - as.numeric(null_ML_st0))
  
  pval_int <- pchisq(abs(LRT_stat_int), df=1, lower.tail = FALSE)
  pval_st1 <- pchisq(abs(LRT_stat_st1), df=1, lower.tail = FALSE)
  pval_st0 <- pchisq(abs(LRT_stat_st0), df=1, lower.tail = FALSE)
  
  #output info
  chr <- unique(data$Chrom)
  pos <- unique(data[data$Window == window, ]$Leading_Position)
  snpid <- unique(data[data$Window == window, ]$Leading_SNP)
  lmaf <- unique(data[data$Window == window, ]$Leading_MAF)
  al0 <- unique(data[data$Window == window, ]$Leading_Allele0)
  al1 <- unique(data[data$Window == window, ]$Leading_Allele1)
  r10 <- freqs$rvals[1]
  r1x <- freqs$rvals[2]
  rx0 <- freqs$rvals[3]
  D10 <- freqs$Dvals[1]
  D1x <- freqs$Dvals[2]
  Dx0 <- freqs$Dvals[3]
  
  #Output results 
  out <- data.frame(chr,window,pos,snpid,lmaf,al0,al1,r10,r1x,rx0,D10,D1x,Dx0,   
                    s_int, maxl_int, s_st1, maxl_st1, s_st0, maxl_st0,
                    pval_int, pval_st1, pval_st0)
  
  colnames(out) <- (c("Chrom","Window","Leading_Position","Leading_SNP","Leading_MAF","Leading_Allele0","Leading_Allele1",
                      "r_10", "r_1x_t", "r_x0_t", "D_10","D_1x_t","D_x0_t",
                      "s_t_int","ML_t_int","s_t_site1","ML_t_site1","s_t_site0","ML_t_site0",
                      "Pval_t_int_lrt", "Pval_t_site1_lrt", "Pval_t_site0_lrt"))
  
  output <- rbind(output, out)
  
  return(output)
}




#Optimization function for interpolation (fecundity)
likelihood.out.f<- function(data,window, pxhat, shrinkage=FALSE){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- as.numeric(unique(unlist(data[data$Window == window,3])))
  expected_haps <- c(0,1,2,3)
  
  #variables
  fem_counts <- NULL
  male_counts <- NULL
  adj_fem_counts <- NULL
  adj_male_counts <- NULL
  freqs0_v <- NULL
  freqs1_v <- NULL
  freqs0_t <- NULL
  freqs1_t <- NULL
  px_hat <- pxhat
  
  #get counts
  window_counts_via <- data[data$Window == window,3:5]
  window_counts_tot <- data[data$Window == window,c(3,10:11)]
  missing_haplotypes_v <- setdiff(expected_haps, window_counts_via$haplotype)
  missing_haplotypes_t <- setdiff(expected_haps, window_counts_tot$haplotype)
  missing_data_v <- data.frame(haplotype = missing_haplotypes_v,
                               female_counts = rep(0, length(missing_haplotypes_v)),
                               male_counts = rep(0, length(missing_haplotypes_v)))
  missing_data_t <- data.frame(haplotype = missing_haplotypes_t,
                               adj_female_counts = rep(0, length(missing_haplotypes_t)),
                               adj_male_counts = rep(0, length(missing_haplotypes_t)))
  combined_data_v <- rbind(window_counts_via, missing_data_v)
  combined_data_t <- rbind(window_counts_tot, missing_data_t)
  
  fem_counts <- as.numeric(unlist(combined_data_v[order(combined_data_v$haplotype), "female_counts"]))
  male_counts <- as.numeric(unlist(combined_data_v[order(combined_data_v$haplotype), "male_counts"]))
  adj_fem_counts <- as.numeric(unlist(combined_data_t[order(combined_data_t$haplotype), "adj_female_counts"]))
  adj_male_counts <- as.numeric(unlist(combined_data_t[order(combined_data_t$haplotype), "adj_male_counts"]))
  
  #get freqs
  freqs_v <- hap_freqs(fem_counts, male_counts, p = px_hat, shrinkage = shrinkage)
  freqs_t <- hap_freqs(adj_fem_counts, adj_male_counts, p = px_hat, shrinkage = shrinkage)
  freqs0_v <- freqs_v$frq_0
  freqs1_v <- freqs_v$frq_1
  freqs0_t <- freqs_t$frq_0
  freqs1_t <- freqs_t$frq_1
  
  
  #Viabiliity optimize (dynamic intervals)
  ml_v <- suppressWarnings(optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                                    p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                                    freqs0 = freqs0_v, freqs1 = freqs1_v))
  
  if (is.nan(ml_v$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = fem_counts, male_counts = male_counts,
                                                                  freqs0 = freqs0_v, freqs1 = freqs1_v))]))
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = fem_counts, male_counts = male_counts,
                                                                    freqs0 = freqs0_v, freqs1 = freqs1_v))]))
    }
    
    ml_v <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                     p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                     freqs0 = freqs0_v, freqs1 = freqs1_v)
  }
  
  
  #site 1 (leftmost) optimize (single)
  scounts1_v <- c(fem_counts[1]+fem_counts[2],
                  male_counts[1]+male_counts[2],
                  fem_counts[3]+fem_counts[4],
                  male_counts[3]+male_counts[4])
  
  ml_st1_v <- optimize(lik.function.single,
                       scounts1_v,
                       interval = c(-0.5,0.5), 
                       maximum = TRUE)
  
  #site 0 (rightmost) optimize (single)
  
  scounts0_v <- c(fem_counts[1]+fem_counts[3],
                  male_counts[1]+male_counts[3],
                  fem_counts[2]+fem_counts[4],
                  male_counts[2]+male_counts[4])
  
  ml_st0_v <- optimize(lik.function.single,
                       scounts0_v,
                       interval = c(-0.5,0.5), 
                       maximum = TRUE)
  
  
  ### get viability results
  s_int_v <- ml_v$maximum
  maxl_int_v <- ml_v$objective
  s_st1_v <- ml_st1_v$maximum
  maxl_st1_v <- ml_st1_v$objective
  s_st0_v <- ml_st0_v$maximum
  maxl_st0_v <- ml_st0_v$objective
  
  
  
  #Total selection optimize (dynamic intervals)
  ml_t <- suppressWarnings(optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                                    p = px_hat, fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                    freqs0 = freqs0_t, freqs1 = freqs1_t))
  
  if (is.nan(ml_t$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                                                  freqs0 = freqs0_t, freqs1 = freqs1_t))]))
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                                                                    freqs0 = freqs0_t, freqs1 = freqs1_t))]))
    }
    
    ml_t <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                     p = px_hat, fem_counts = adj_fem_counts, male_counts = adj_male_counts,
                     freqs0 = freqs0_t, freqs1 = freqs1_t)
  }
  
  
  #site 1 (leftmost) optimize (single)
  scounts1_t <- c(adj_fem_counts[1]+adj_fem_counts[2],
                  adj_male_counts[1]+adj_male_counts[2],
                  adj_fem_counts[3]+adj_fem_counts[4],
                  adj_male_counts[3]+adj_male_counts[4])
  
  ml_st1_t <- optimize(lik.function.single,
                       scounts1_t,
                       interval = c(-0.5,0.5), 
                       maximum = TRUE)
  
  #site 0 (rightmost) optimize (single)
  
  scounts0_t <- c(adj_fem_counts[1]+adj_fem_counts[3],
                  adj_male_counts[1]+adj_male_counts[3],
                  adj_fem_counts[2]+adj_fem_counts[4],
                  adj_male_counts[2]+adj_male_counts[4])
  
  ml_st0_t <- optimize(lik.function.single,
                       scounts0_t,
                       interval = c(-0.5,0.5), 
                       maximum = TRUE)
  
  
  ### get total selection results
  s_int_t <- ml_t$maximum
  maxl_int_t <- ml_t$objective
  s_st1_t <- ml_st1_t$maximum
  maxl_st1_t <- ml_st1_t$objective
  s_st0_t <- ml_st0_t$maximum
  maxl_st0_t <- ml_st0_t$objective
  
  
  
  ### Calculate fecundity selection coefficients
  p<-pxhat
  q<-1-pxhat
  s_int_f <- (s_int_t-s_int_v)/((1+q*s_int_v)*(1-p*s_int_v))
  s_st1_f <- (s_st1_t-s_st1_v)/((1+q*s_st1_v)*(1-p*s_st1_v))
  s_st0_f <- (s_st0_t-s_st0_v)/((1+q*s_st0_v)*(1-p*s_st0_v))
  
  #output info
  chr <- unique(data$Chrom)
  pos <- unique(data[data$Window == window, ]$Leading_Position)
  snpid <- unique(data[data$Window == window, ]$Leading_SNP)
  lmaf <- unique(data[data$Window == window, ]$Leading_MAF)
  al0 <- unique(data[data$Window == window, ]$Leading_Allele0)
  al1 <- unique(data[data$Window == window, ]$Leading_Allele1)

  #Output results 
  out <- data.frame(chr,window,pos,snpid,lmaf,al0,al1,
                    s_int_f,s_st1_f, s_st0_f)
  
  colnames(out) <- (c("Chrom","Window","Leading_Position","Leading_SNP","Leading_MAF","Leading_Allele0","Leading_Allele1",
                      "s_f_int","s_f_site1","s_f_site0"))
  
  output <- rbind(output, out)
  
  return(output)
}