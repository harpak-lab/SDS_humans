## 3. SAMPLE-LEVEL FILTERING OF UK BIOBANK
## updated 5/1/24
## Jared M. Cole
##
## Extract sample IDs that meet QC criteria  
###########################################################################

setwd("/Documents/Projects/SAS/humans/UKB/metadata/")
library(data.table)

# import sample ID list from UKB BGEN files
samples <- read.table("ukb61666_hap_chr1_v2_s487220.sample", 
	sep=" ", header = TRUE, check.names=FALSE)
samples <- samples[-1,-c(2,3)]
names(samples)[1] <- "eid"

# load extracted metadata for relevant fields
md <- read.table("UKB_metadata_subset_1-24.csv", sep=",", 
	header = TRUE, check.names=FALSE)

### QC filtering

# remove heterozygote and missingness outliers (22027, 1=outlier)
md$`22027-0.0`[is.na(md$`22027-0.0`)] <- 0 #set NAs to 0
md.ft.1 <- md[md$`22027-0.0` == 0,]
rm(md)

# remove sex chromosome aneuploids (22019, 1=aneuploids)
md.ft.1$`22019-0.0`[is.na(md.ft.1$`22019-0.0`)] <- 0
md.ft.2 <- md.ft.1[md.ft.1$`22019-0.0` == 0,]
rm(md.ft.1)

# remove sex mismatches (21 and 22001, female=0, male =1)
md.ft.3 <- md.ft.2[!is.na(md.ft.2$`31-0.0`),] #remove samples with NAs for either column
md.ft.3 <- md.ft.3[!is.na(md.ft.3$`22001-0.0`),]
md.ft.3 <- md.ft.3[md.ft.3$`31-0.0` == md.ft.3$`22001-0.0`,]
rm(md.ft.2)

# extract only samples used in PCA, controlled for relatedness (field 22020, 1=yes)
md.ft.4 <- md.ft.3[!is.na(md.ft.3$`22020-0.0`),]
rm(md.ft.3)

# extract only white, British subsample (field, 22006, 1=Caucasian)
md.ft.5 <- md.ft.4[!is.na(md.ft.4$`22006-0.0`),]
rm(md.ft.4)


#### Filtering for lifetime reproductive success (LRS)

#split by sex
m <- subset(md.ft.5, `31-0.0` == 1)
f <- subset(md.ft.5, `31-0.0` == 0)

#maximum number of offspring across assessments (2405 for males, 2734 for females)
m$LRS <- apply(m[,c("2405-0.0","2405-1.0","2405-2.0","2405-3.0")],1,
	function(x)max(x,na.rm=T))
f$LRS <- apply(f[,c("2734-0.0","2734-1.0","2734-2.0","2734-3.0")],1,
	function(x)max(x,na.rm=T))

#removing -3 (didn't answer), -1 (do not know), and -Inf (not recorded, NA)
m <- m[!m$LRS == -Inf,]
f <- f[!f$LRS == -Inf,]
m <- m[m$LRS >= 0 ,]
f <- f[f$LRS >= 0 ,]

#check which samples report more children at earlier assessment points than later ones and remove
m$x <- ifelse(m$`2405-3.0`< m$`2405-2.0` | m$`2405-3.0`<m$`2405-1.0` | m$`2405-2.0`<m$`2405-1.0` | 
	m$`2405-2.0`<m$`2405-0.0` | m$`2405-3.0`<m$`2405-0.0`| m$`2405-1.0`<m$`2405-0.0`, "yes", "no")
f$x <- ifelse(f$`2734-3.0`< f$`2734-2.0` | f$`2734-3.0`<f$`2734-1.0` | f$`2734-2.0`<f$`2734-1.0` | 
	f$`2734-2.0`<f$`2734-0.0` | f$`2734-3.0`<f$`2734-0.0`| f$`2734-1.0`<f$`2734-0.0`, "yes", "no")

m.1 <- m[m$x == "no" | is.na(m$x),]
f.1 <- f[f$x == "no" | is.na(f$x),]

rm(md.ft.5)
rm(m)
rm(f)

#combine samples
md.ft.6 <- rbind(m.1,f.1)

rm(m.1)
rm(f.1)

#remove individuals under 45 (field 21022) and with more than 20 offspring
md.ft.7 <- md.ft.6[md.ft.6$LRS < 20 &  md.ft.6$`21022-0.0` >= 45, ]
rm(md.ft.6)

#remove withdraws from the dataset (since 5/2023)
wd <- read.table("withdraws_2023.txt", sep="", header = FALSE, 
	check.names=FALSE)
names(wd)[1] <- "eid"
md.ft.8 <- md.ft.7[!md.ft.7$eid %in% wd$eid,]
rm(md.ft.7)

# match samples IDs from metadata and UKB BGEN
md.ft.9 <- md.ft.8[md.ft.8$eid %in% samples$eid,]

# export filtered metadata
write.table(md.ft.9, "UKB_metadata_filtered_1-24.txt", row.names = FALSE, 
	sep="\t", quote = FALSE)

#export samples list for PLINK and ML pipeline
ex_samples <- samples[samples$eid %in% md.ft.9$eid,]
ex_samples.1<-merge(md.ft.9, ex_samples, by="eid")
ex_samples.2 <- ex_samples.1[,c(1,22,20)]
ex_samples.2$ID_2 <- ex_samples.2$eid
ex_samples.2$missing <- 0
names(ex_samples.2)[3] <- "children"
ex_samples.3 <- ex_samples.2[,c(1,3)]
names(ex_samples.2)[1] <- "ID_1"
ex_samples.4 <- ex_samples.2[,c(1,4)]

#xport
write.table(ex_samples.3, "UKB_sample_LRS_1-24.txt", row.names = FALSE, 
	sep="\t", quote = FALSE)
write.table(ex_samples.4, "UKB_filtered_samples_list_1-24.txt", 
	row.names = FALSE, sep="\t",
            quote = FALSE, col.names = FALSE)
