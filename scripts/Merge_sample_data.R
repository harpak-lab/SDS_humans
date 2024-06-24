#### Merge samples with LRS

#Remove scientific notation
options(scipen = 999)
options(warn=1)

#Argument, libraries, and settings
suppressMessages(library("optparse")) #load opt parser library
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="input file", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="samples", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", 
    call.=FALSE)
}
if (is.null(opt$samples)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (samples file).n", 
    call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", 
    call.=FALSE)
}

#Load LRS file
LRS <- read.table(opt$file, sep="", header = TRUE, 
  check.names=FALSE)
LRS <- rbind(c(0,"C"), LRS)
names(LRS) <- c("ID_1","children")

#Load samples
samples <- read.table(opt$samples, sep="", header = TRUE, 
  check.names=FALSE, stringsAsFactors = FALSE)

#Merge
merged_samples <- merge(samples, LRS, by = "ID_1", sort = F)

#export
write.table(merged_samples, opt$output, row.names = FALSE, sep=" ", 
  quote = FALSE)