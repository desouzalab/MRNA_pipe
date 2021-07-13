#======================================================================================================
# BackSpin Results in R
#======================================================================================================
suppressMessages(library(argparse))
suppressMessages(library(scater))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(hydroGOF))
suppressMessages(library(glmnet))
suppressMessages(library(StatMatch))
suppressMessages(library(Rtsne))
suppressMessages(library(fpc))
suppressMessages(library(GA))
suppressMessages(library(MASS))
suppressMessages(library(session))
suppressMessages(library(Matrix))
suppressMessages(library(vegan))
suppressMessages(library(data.table))
suppressMessages(library(reshape))
suppressMessages(library(abind))
suppressMessages(library(drc))
#======================

RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing what")
parser$add_argument("--data_output_directory", type="character", help="Path to the cluster data output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE, recursive=TRUE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE, recursive=TRUE)

all_cef_result_files <- list.files(args$input_directory, pattern="*.cef")
print(all_cef_result_files)

# Results in R
for (c in 1:length(all_cef_result_files)){
  print(c)
  print(all_cef_result_files)
  ### Read cef document
  cef=read.csv(file.path(args$input_directory, all_cef_result_files[c]), header=F, sep="\t")
  cef=as.matrix(cef)
  # # ordered count matrix: cef.dat
  cef.dat=cef[11:nrow(cef), 10:ncol(cef)]
  cef.dat=t(apply(cef.dat, 1, as.numeric))
  colnames(cef.dat)=as.vector(cef[2, 10:ncol(cef)])
  rownames(cef.dat)=as.vector(cef[11:nrow(cef), 1])
  # results for gene cluster
  cef.gene.k=cef[11:nrow(cef), 3:8]
  cef.gene.k=t(apply(cef.gene.k, 1, as.numeric))
  rownames(cef.gene.k)=as.vector(cef[11:nrow(cef), 1])
  colnames(cef.gene.k)=1:6
  #results for cell cluster
  cef.cell.k=t(cef[4:9, 10:ncol(cef)])
  cef.cell.k=t(apply(cef.cell.k, 1, as.numeric))
  rownames(cef.cell.k)=as.vector(cef[2, 10:ncol(cef)])
  colnames(cef.cell.k)=paste("CellType_L",1:6,sep="") # cluster results under iteration 1:6
  cef.cell.6 = cef.cell.k[,6]
  #write.csv(cef.cell.k, file=paste0(data_outdir,"/clustersBackSpin_k",c,"_",args$name_dataset,".csv"))
  write.csv(cef.cell.6, file=paste0(data_outdir,"/clustersBackSpin6_",c,"_",args$name_dataset,".csv"))
  #write.csv(cef.gene.k, file=paste0(data_outdir,"/clustersBackSpin_forGene_",c,"_",args$name_dataset,".csv"))
}
print("DONE")