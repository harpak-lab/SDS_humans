## Run ABC
# updated 6/11/24
# author: Carly B Scott

##############################################

#ABC procedure 
library(abc)
library(ggExtra)
library(tidyverse)

abc_res = read.table("FinalABC_50kSims.output", header = F) # the COMBINED output of your simulations
colnames(abc_res) <- c("S_prior", "F_prior", "SSE")

#plot priors, check
par(mfrow=c(1,2))
hist(log10(abc_res$S_prior), main = "S Prior", xlab = "log10(s)")
hist(log10(abc_res$F_prior), main = "F Prior", xlab = "log10(F)")
dev.off()


#Run ABC on log transformed values (doesn't matter, but makes viz easier)
abc_res[,1] <- log10(abc_res[,1])
abc_res[,2] <- log10(abc_res[,2])
#aim for a target of zero - ideally there would be zero error between observed and simulated values. 
testabc <- abc(target = 0, param = abc_res[,1:2], sumstat = abc_res[,3], tol = .01, method = "rejection")
summary(testabc) #check the summary output
abc.df <- data.frame(testabc$unadj.values)
colnames(abc.df) = c("accepted_s", "accepted_F")

#find the maximum density of the joint distribution by doing a kernel density estimation
smooth <- kde2d(x=abc.df$accepted_s, y=abc.df$accepted_F, n = 1000)
which.max(smooth$z)
getrowcol <- 823412/1000
getrowcol
xmode = smooth$x[412] #this comes from the which.max value. We know the dimensions of our raster, so can find location in this way.
#I did this manually, double check it below for your data. 
ymode= smooth$y[823]
image(smooth)
points(xmode, ymode, pch = 19, cex = 2)

write.table(abc.df, "ABC_Results_Top1Perc.tsv", sep = "\t", col.names = T, row.names = F)

########################################################
##### CHECK SENSITIVITY TO THRESHHOLDS BELOW###########
#This will overwrite some of the objects you've already created above. 

#downsample 
#get max density 
outdf <- data.frame()
for(i in rep(seq(.05, 1, .1), 50)){
  abc_res_dwn <- abc_res[sample(rownames(abc_res), i*nrow(abc_res), replace = T),]
  testabc <- abc(target = 0, param = abc_res_dwn[,1:2], sumstat = abc_res_dwn[,3], tol = .01, method = "rejection")
  abc.df <- data.frame(testabc$unadj.values)
  colnames(abc.df) = c("accepted_s", "accepted_F")
  smooth <- kde2d(x=abc.df$accepted_s, y=abc.df$accepted_F, n = 100)
  n <- which.max(smooth$z)
  getrowcol <- n/100
  getrowcol
  xmode = smooth$x[(getrowcol-floor(getrowcol))*100]
  smean=mean(abc.df$accepted_s)
  fmean=mean(abc.df$accepted_F)
  ymode= smooth$y[floor(getrowcol)]
  outdf <- rbind(outdf, data.frame(xmode=xmode,ymode=ymode,smean=smean, fmean=fmean,dwnsmp=i))
}
colnames(outdf) <- c("S_posterior_mode", "F_posterior_mode", "S_posterior_mean", "F_posterior_mean","DownsamplePercent")

write.table(outdf, "DownSampling_vs_ABCestimates.txt", row.names = F, quote = F)

outdf <- data.frame()
for(i in seq(.005, .2, .001)){
  #abc_res_dwn <- abc_res[sample(rownames(abc_res), i*nrow(abc_res), replace = T),]
  testabc <- abc(target = 0, param = abc_res_dwn[,1:2], sumstat = abc_res_dwn[,3], tol = i, method = "rejection")
  abc.df <- data.frame(testabc$unadj.values)
  colnames(abc.df) = c("accepted_s", "accepted_F")
  smooth <- kde2d(x=abc.df$accepted_s, y=abc.df$accepted_F, n = 100)
  n <- which.max(smooth$z)
  getrowcol <- n/100
  getrowcol
  xmode = smooth$x[(getrowcol-floor(getrowcol))*100]
  smean=mean(abc.df$accepted_s)
  fmean=mean(abc.df$accepted_F)
  ymode= smooth$y[floor(getrowcol)]
  outdf <- rbind(outdf, data.frame(xmode=xmode,ymode=ymode,smean=smean, fmean=fmean,dwnsmp=i))
}
colnames(outdf) <- c("S_posterior_mode", "F_posterior_mode", "S_posterior_mean", "F_posterior_mean","Acceptance_Threshold")
write.table(outdf, "Threshold_vs_ABCestimates.txt", row.names = F, quote = F)
