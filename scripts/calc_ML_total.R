# Calculate ML for haplotype data
# SAS total selection

# updated - 1/10/24

#Remove scientific notation
#options(scipen = 999)  
options(warn=1)

#Argument, libraries, and settings
suppressMessages(library("optparse")) #load opt parser library
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="haplotype counts file", metavar="character"),
  make_option(c("-s", "--snpfile"), type="character", default=NULL, 
              help="snp file", metavar="character"),
  make_option(c("-w", "--windowsize"), type="character", default=5, 
              help="window size [default= %default]", metavar="numeric"),
  make_option(c("-p", "--pxhat"), type="character", default=0.2, 
              help="assumed minor allele frequency at site x  [default= %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input haplotype file).n", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", call.=FALSE)
}

### Import functions
source('/scratch/ukb/scripts/ML_functions.R')

#input haplotype counts and count matrices
haplotype_counts <- fread(opt$file, sep="\t", header = TRUE, check.names=FALSE)
snps <- fread(opt$snpfile, sep="\t", header = TRUE, check.names=FALSE)

#Window size
win <- as.numeric(opt$windowsize)

#px_hat
pxhat <- as.numeric(opt$pxhat)   

#Likelihood calculations
winlist <- as.character(unique(haplotype_counts$Window))
ml_result <- lapply(winlist, function(window) likelihood.out.t(haplotype_counts, window, 
                                                               pxhat = pxhat, shrinkage=FALSE))
ml_output <- do.call(rbind, ml_result)
ml_output$Window <- as.numeric(ml_output$Window)
ml_output[order(ml_output$Window),] -> ml_output

# Get SNP info for windows
snp_mapping <- data.frame(Site1_SNP = snps$SNP[-nrow(snps)], 
                          Site0_SNP = snps$SNP[-1])
ml_output1 <- merge(ml_output, snp_mapping, by.x = "Leading_SNP", by.y = "Site1_SNP", all.x = TRUE)
ml_output1 <- merge(ml_output1 , snps, by.x = "Site0_SNP", by.y = "SNP", all.x = TRUE, suffixes = c("", "_Site0"))
ml_output1 <- within(ml_output1, {
  Site1_Position <- Leading_Position
  Site1_SNP <- Leading_SNP
  Site1_MAF <- Leading_MAF
  Site1_Allele0 <- Leading_Allele0
  Site1_Allele1 <- Leading_Allele1
  Site0_Position <- Locus
  Site0_MAF <- MAF
  Site0_Allele0 <- Allele0
  Site0_Allele1 <- Allele1
})

ml_output1 <- ml_output1[, c("Chrom", "Window", "Site1_Position", "Site1_SNP", "Site1_MAF", "Site1_Allele0", "Site1_Allele1",
                                 "Site0_Position", "Site0_SNP", "Site0_MAF", "Site0_Allele0", "Site0_Allele1", 
                              "r_10", "r_1x_t", "r_x0_t","D_10", "D_1x_t","D_x0_t","s_t_int","ML_t_int","s_t_site1","ML_t_site1", "s_t_site0", "ML_t_site0",
                              "Pval_t_int_lrt",   "Pval_t_site1_lrt", "Pval_t_site0_lrt")]
ml_output1$Window <- as.numeric(ml_output1$Window)
ml_output1[order(ml_output1$Window),] -> ml_output1

#Output
write.table(ml_output1, paste0(opt$output, ".t.maxl"), row.names = FALSE, sep="\t")
