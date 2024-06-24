# Get Final results from bootstrapping and ML
# SAS viability

# updated - 2/8/24

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
    all.=FALSE)
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
out1 <- lapply(snpids,get_boot_p.v,bootdist=boot,est=ml_output)
p_boot_dist <- do.call(rbind, out1)


#combine with ML results
f_output <- merge(ml_output, p_boot_dist, by="Site1_SNP")
f_output <- f_output %>% select(Chrom,Window,Site1_Position,
  Site1_SNP,Site1_MAF,Site1_Allele0,Site1_Allele1,
  Site0_Position,Site0_SNP,Site0_MAF,Site0_Allele0,
  Site0_Allele1,r_10,r_1x,r_x0,D_10,D_1x,D_x0,
  s_int_viability = s_v_int, ML_int_viability = ML_v_int, SE_v, Z_v, 
  Pval_int_lrt_viability = Pval_v_int_lrt, Pval_int_boot_viability = P_v,
  s_site1_viability =  s_v_site1, ML_site1_viability = ML_v_site1, 
  Pval_site1_lrt_viability = Pval_v_site1_lrt,
  s_site0_viability =  s_v_site0, ML_site0_viability = ML_v_site0, 
  Pval_site0_lrt_viability = Pval_v_site0_lrt, Outside_boot_via=Outside)
f_output[order(f_output$Window),] -> f_output

#Output
write.table(f_output, paste0(opt$output, ".v.result"), 
  row.names = FALSE, sep="\t")
