# Bootstrapping for UKB haplotype data
# SAS total 

# updated - 3/15/24

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
  make_option(c("-m", "--male"), type="character", default=NULL, 
              help="male hap matrix", metavar="character"),
  make_option(c("-g", "--female"), type="character", default=NULL, 
              help="female hap matrix", metavar="character"),
  make_option(c("-w", "--windowsize"), type="character", default=2, 
              help="window size [default= %default]", metavar="numeric"),
  make_option(c("-b", "--bootreps"), type="character", default=100, 
              help="bootstrap samples [default= %default]", metavar="numeric"),
  make_option(c("-p", "--pxhat"), type="character", default=0.2, 
              help="assumed minor allele frequency at site x  [default= %default]", 
              metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input haplotype file).n", 
    call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", 
    call.=FALSE)
}

### Import functions
source('/scratch/ukb/scripts/ML_functions.R')

#input haplotype counts
haplotype_counts <- fread(opt$file, sep="\t", header = TRUE, 
  check.names=FALSE)
male_m <- readRDS(file = opt$male)
female_m <- readRDS(file = opt$female)

#Window size
win <- as.numeric(opt$windowsize)

#px_hat
pxhat <- as.numeric(opt$pxhat)   

# Bootstrap
boots = vector("list", length = as.numeric(opt$bootreps))

for (i in 1:as.numeric(opt$bootreps)) {
  bt_counts <- bootstrap.t(haplotype_counts, 
                           m = male_m, f = female_m, pxhat = pxhat)
  bt_counts$bootrep <- i
  bt_counts[order(bt_counts$Window),] -> bt_counts
  boots[[i]] <- bt_counts 
}
boots = do.call(rbind, boots)
boots$Window <- as.numeric(as.character(boots$Window))
boots <- boots[order(boots$bootrep, boots$Window), ]


#output
write.table(boots, paste0(opt$output, ".t.boot"), sep="\t", 
            row.names = FALSE, quote = FALSE)