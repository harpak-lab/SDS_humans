#Transpose haplotype data, add sample info
## 5/11/23, Jared Cole
##########################################

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
  make_option(c("-s", "--samples"), type="character", default=NULL, 
              help="samples", metavar="character"),
  make_option(c("-v", "--variantfile"), type="character", default=NULL, 
              help="variant file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input haplotype file).n", call.=FALSE)
}
if (is.null(opt$samples)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (samples file).n", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", call.=FALSE)
}
if (is.null(opt$variantfile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (variant file).n", call.=FALSE)
}



#Load haplotype block PLINK files (add in filename from opts)
haploblock <- fread(opt$file, sep=" ", header = FALSE, check.names=FALSE)

#Load variant data and get SNP list (variant file is from snps_file.sh output)
variants <- read.table(opt$variantfile, sep="", header = TRUE, check.names=FALSE)
analysis_vars <- as.data.frame(haploblock[,c(2,4,5)])
names(analysis_vars) <- c("SNP","Allele0","Allele1")
snps <- inner_join(variants, analysis_vars, by="SNP") %>%
  select("SNP", "Locus"="POS", "MAF", "Allele0","Allele1") %>% unique()

#Load haplotype and sample
samples <- read.table(opt$samples, sep="", header = TRUE, check.names=FALSE, stringsAsFactors = FALSE)
samples <- samples[-1,] 
samples %>% select("sample" = "ID_2", "sex", "children") -> samples

#Transpose dataframe
hap_test <- as.data.frame(t(haploblock))
hap_test <- hap_test[-c(1,2,4,5),]
colnames(hap_test) <- as.character(unlist(hap_test[1,]))
hap_test = hap_test[-1, ]

#Add columns with samples and sex 
sample_list <- c(rep(samples$sample, each=2))
sex_list <- c(rep(samples$sex, each=2))
lrs_list <- c(rep(samples$children, each=2))
hap_test %>%
  mutate(sample = sample_list, sex = sex_list, lrs = lrs_list) %>%
  select (sample, sex, lrs, everything()) %>%
  mutate_all(as.character) -> hap_test

#export
write.table(hap_test, paste0(opt$output, ".hap"), row.names = FALSE, sep="\t")
write.table(snps, paste0(opt$output, ".snps"), row.names = FALSE, sep="\t")
