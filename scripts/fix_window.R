# Fix window column for final output

#Remove scientific notation
options(scipen = 999)  
options(warn=1)

#Argument, libraries, and settings
library("optparse") #load opt parser library
library("dplyr")

option_list = list(
  make_option(c("-f", "--file1"), type="character", default=NULL, 
              help="Concatentated Result file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Result output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file1)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", 
  	call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file name).n", 
  	call.=FALSE)
}

#results
result <- read.table(paste0(opt$file1), sep="", header = TRUE, 
	check.names=FALSE)

#remove duplicated SNPs
result <- result[!duplicated(result[, 3]), ]

#windows
winds <- nrow(result)

#order
result<-result[order(result$Site1_Position),]  

#put correct windows
result$Window <- rep(1:winds)

#out
write.table(result, paste0(opt$output, ".combined"), 
	row.names = FALSE, sep="\t")