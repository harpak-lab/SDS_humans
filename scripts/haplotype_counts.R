## Calculate Haplotype Counts table for UKB haplotype data
## 1/9/24, Jared Cole
##########################################################

#Remove scientific notation
options(scipen = 999)  
options(warn=1)

#Argument, libraries, and settings
suppressMessages(library("optparse")) #load opt parser library
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="haplotype file", metavar="character"),
  make_option(c("-w", "--windowsize"), type="character", default=5, 
              help="window size [default= %default]", metavar="numeric"),
  make_option(c("-c", "--chrom"), type="character", 
              help="chromosome", metavar="numeric"),
  make_option(c("-v", "--snps"), type="character", default=NULL, 
              help="snp file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input haplotype file).n", call.=FALSE)
}
if (is.null(opt$snps)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (snp file).n", call.=FALSE)
}
if (is.null(opt$chrom)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (variant file).n", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", call.=FALSE)
}

### Import functions
source('/scratch/ukb/scripts/ML_functions.R')

#Load haplotype block PLINK files (add in filename from opts)
hap_test <- fread(opt$file, sep="\t", header = TRUE, check.names=FALSE)
hap_test <- as.data.frame(hap_test) #to get hap counts as dataframe

#chromosome
chr <- as.numeric(opt$chrom)

#Load variant data and get SNP list (variant file is from snps_file.sh output)
snps <- read.table(opt$snps, sep="", header = TRUE, check.names=FALSE, stringsAsFactors = FALSE)

### Windowed analysis
cols <- 4:ncol(hap_test)
win <- as.numeric(opt$windowsize)   #Window size
wins <- getWindows(cols,win,(win-1))

# Make sure each window has the same length
wins <- wins[lengths(wins)==win]

# Get haplotypes for males and females
males <- lapply(wins,gethaps,input=hap_test,s=1) #1 for males
females <- lapply(wins,gethaps,input=hap_test,s=2) #2 for females

#Remove Ambiguous Sites (with ?)
#males <- lapply(males, function(x) x[!grepl("?", x)])
#females <- lapply(females, function(x) x[!grepl("?", x)])

# Get LRS for male and female haplotypes
males_lrs <- lapply(wins,getlrs,input=hap_test,s=1) #1 for males
females_lrs <- lapply(wins,getlrs,input=hap_test,s=2) #2 for females

# Get frequency matrices for fecundity and total selection
male_matrix <- haplrs_matrix(males,males_lrs,1)
female_matrix <- haplrs_matrix(females,females_lrs,1)

# Validate haplotype columns across two matrices
male_val <- haplrs_validate(male_matrix,female_matrix)
female_val <- haplrs_validate(female_matrix,male_matrix)

#Count the haplotypes and convert to base 10
males_counts <- haplocounts(males, males_lrs, 1) 
females_counts <- haplocounts(females, females_lrs, 2)

hapcounts <- data.frame(Window = 1:length(wins))
colnamevec <- c("haplotype", "count", "lrs")
hapcounts[ ,colnamevec] <- NA

for (i in 1:length(females_counts)) {
  females_counts[[i]]$Window <- as.numeric(names(females_counts[i]))
}
for (i in 1:length(males_counts)) {
  males_counts[[i]]$Window <- as.numeric(names(males_counts[i]))
}

female_haps <- do.call(rbind, females_counts)
male_haps <- do.call(rbind, males_counts)
female_haps$sex <- 2
male_haps$sex <-1

#Create male and female dataframe for haplotype counts
merge(hapcounts, female_haps, by=c("Window","haplotype"), all=TRUE) -> merged_df_f
merge(hapcounts, male_haps,  by=c("Window","haplotype"), all=TRUE) -> merged_df_m
merged_df_f %>%
  select(-count.x,-lrs.x) -> merged_df_f
merged_df_m %>%
  select(-count.x,-lrs.x) -> merged_df_m
na.omit(merged_df_f) -> merged_df_f
na.omit(merged_df_m) -> merged_df_m

#Merge male and female dataframes
merge(merged_df_f, merged_df_m, by=c("Window","haplotype"), all = TRUE) -> merged_df
haplotype_counts <- merged_df %>% 
  select(Window, haplotype, female_counts = count.y.x, male_counts = count.y.y, female_lrs = lrs.y.x,
         male_lrs = lrs.y.y)
haplotype_counts[is.na(haplotype_counts)] = 0

#fecundity calculations 
haplotype_counts$female_lrs2 <- haplotype_counts$female_lrs/2 #halve the total LRS
haplotype_counts$male_lrs2 <- haplotype_counts$male_lrs/2
haplotype_counts <- haplotype_counts %>% group_by(Window) %>% #get avg LRS for both sexes
  mutate(aF = sum(female_lrs2)/sum(female_counts))
haplotype_counts <- haplotype_counts %>% group_by(Window) %>% 
  mutate(aM = sum(male_lrs2)/sum(male_counts))
haplotype_counts$adj_female_counts <- haplotype_counts$female_lrs/(2 * haplotype_counts$aF) #weight haplotype
haplotype_counts$adj_male_counts <- haplotype_counts$male_lrs/(2 * haplotype_counts$aM)

#Get SNP ids and MAFs for leading SNP in window
haplotype_counts$Chrom <- chr
haplotype_counts$Locus <- sapply(haplotype_counts$Window, getsnppos, leading=1, input=hap_test)
haplotype_counts <- merge(snps, haplotype_counts, by="Locus")
haplotype_counts <- haplotype_counts %>%
  select(Chrom, Window, haplotype, female_counts, male_counts, female_lrs, male_lrs, female_lrs2, male_lrs2, adj_female_counts,
         adj_male_counts, Leading_Position = Locus, Leading_SNP = SNP, Leading_MAF = MAF, Leading_Allele0 = Allele0, Leading_Allele1 = Allele1)

#Output
write.table(haplotype_counts, paste0(opt$output, ".counts"), row.names = FALSE, sep="\t")
saveRDS(male_val, file = paste0(opt$output, ".malematrix.rds"))
saveRDS(female_val, file = paste0(opt$output, ".femalematrix.rds"))