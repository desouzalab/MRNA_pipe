#======================================================================================================
# cef transformation for preprocessed data as input of BackSpin
#======================================================================================================

#======================

RSCRIPT <- "Rscript"

#======================
# arguments
#======================
suppressMessages(library(argparse))


RSCRIPT <- "Rscript"
# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--data_output_directory", type="character", help="Path to the cef ssRNA data output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE, recursive=TRUE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$input_directory, pattern="*.csv*")
print(all_preprocessed_ssRNASeq_files)


for (c in 1:length(all_preprocessed_ssRNASeq_files)){
  print(c)
  print(all_preprocessed_ssRNASeq_files)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_preprocessed_ssRNASeq_files[c]),row.names = 1)
  data=na.omit(data)
  print("  ...read")
  #transforming to cef document
  cef.dat=data
  output.cef=paste(data_outdir,"/ceftype","_",args$name_dataset,".cef", sep="")
  cef=cef.dat
  cef=rbind(gene="", data.frame(well="", cef) )
  cef.head=c("CEF",  "0", "1","1", nrow(cef.dat), ncol(cef.dat), "0")
  write.table(matrix(cef.head, nrow=1),output.cef, sep="\t", row.names=F, col.names=F,quote=F)
  write.table(cef, output.cef, sep="\t", row.names=T, col.names=NA,quote=F, append = T)
  rm(list = ls())

}
print("DONE")
#======================================================================================================