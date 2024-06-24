# Get Final results from bootstrapping and ML
# SAS total

# updated - 3/28/24

#Remove scientific notation
#options(scipen = 999)  
options(warn=1)

#Argument, libraries, and settings
suppressMessages(library("optparse")) #load opt parser library
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))

option_list = list(
  make_option(c("-m", "--mlfile"), type="character", default=NULL, 
              help="ML results file (.maxl)", metavar="character"),
  make_option(c("-b", "--bootfile"), type="character", default=NULL, 
              help="bootstrap file (.boot)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$mlfile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input maxl file).n", 
    call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", 
    call.=FALSE)
}

### Import functions
source('/scratch/ukb/scripts/ML_functions.R')

#ML file
ml_output <- fread(opt$mlfile, sep="\t", header = TRUE, 
  check.names=FALSE)

## UPLOAD BOOTS FILE
boot <- fread(opt$bootfile, sep="\t", header = TRUE, 
  check.names=FALSE)

#get p-values from bootstrapping
snpids <- as.vector(unique(ml_output$Site1_SNP))
p_boot_dist <- data.frame()
out1 <- lapply(snpids,get_boot_p.t,bootdist=boot,est=ml_output)
p_boot_dist <- do.call(rbind, out1)


#combine with ML results
f_output <- merge(ml_output, p_boot_dist, by="Site1_SNP")
f_output <- f_output %>% select(Chrom,Window,Site1_Position,
  Site1_SNP,Site1_MAF,Site1_Allele0,Site1_Allele1,
  Site0_Position,Site0_SNP,Site0_MAF,Site0_Allele0,
  Site0_Allele1,r_10,r_1x_t,r_x0_t,D_10,D_1x_t,D_x0_t,
  s_int_total = s_t_int, ML_int_total = ML_t_int, SE_t, 
  Z_t, Pval_int_lrt_total = Pval_t_int_lrt, Pval_int_boot_total = P_t,
  s_site1_total =  s_t_site1, ML_site1_total = ML_t_site1, 
  Pval_site1_lrt_total = Pval_t_site1_lrt,
  s_site0_total =  s_t_site0, ML_site0_total = ML_t_site0, 
  Pval_site0_lrt_total = Pval_t_site0_lrt, Outside_boot_via=Outside)
f_output[order(f_output$Window),] -> f_output

#Output
write.table(f_output, paste0(opt$output, ".t.result"), row.names = FALSE, sep="\t")
