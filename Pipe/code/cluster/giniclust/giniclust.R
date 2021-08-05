
#======================================================================================================
#======================================================================================================
# CLUSTER - GiniClust
#======================================================================================================

#======================
# argparse library
#======================
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

RSCRIPT <- "Rscript"


#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--data_output_directory", type="character", help="Path to the cluster data output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

parser$add_argument("--input_raw", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE, recursive=TRUE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$input_directory, pattern="*.csv*")
print(all_preprocessed_ssRNASeq_files)

raw_data <- list.files(args$input_raw, pattern="*.csv*")
print(raw_data)


for (c in 1:length(all_preprocessed_ssRNASeq_files)){
  print(c)
  print(all_preprocessed_ssRNASeq_files)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_preprocessed_ssRNASeq_files[c]),row.names = 1)
  data=na.omit(data)
  ExprM.RawCounts=na.omit(read.csv(file.path(args$input_raw, raw_data[c]),row.names = 1))
  print("  ...read")
  #parameters
  minCellNum           = 3                                                # filtering, for at least expressed in how many cells
  minGeneNum           = 2000                                             # filtering, for at least expressed in how many genes
  expressed_cutoff     = 1                                                # filtering, for raw counts
  log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression   
  log2.expr.cutoffh    = 20                                               # cutoff for range of gene expression 
  Gini.pvalue_cutoff   = 0.0001                                           # fiting, Pvalue, control how many gene finally used.
  Norm.Gini.cutoff     = 1                                                # fiting, NomGini, control how many gene finally used, 1 means not used.
  span                 = 0.9                                              # parameter for LOESS fitting
  outlier_remove       = 0.75                                             # parameter for LOESS fitting
  Gamma                = 0.9                                              # parameter for clustering
  diff.cutoff          = 1                                                # MAST analysis, filter gene don't have high log2_foldchange to reduce gene num
  lr.p_value_cutoff    = 1e-5                                             # MAST analysis, pvalue cutoff to identify differential expressed gene
  CountsForNormalized  = 100000                                           
  rare_p               = 0.05                                             # propostion of cell number < this value will be considered as rare cell clusters.
  perplexity           = 30
  eps                  = 0.5                                              # parameter for DBSCAN
  MinPts               = 3                                                # parameter for DBSCAN
  ExprM.RawCounts.filter <- data
  

  source("~/projects/def-cdesouza/Lab/GiniClust/Rfunction/GiniClust_Fitting.R")#Upload the function files and give the right path to call them
  GeneList.final = GiniClust_Fitting(data.type = 'RNA-seq',ExprM.RawCounts.filter=data,out.folder=data_outdir,exprimentID=c)
  print("done fitting")
  #clustering
  source("~/projects/def-cdesouza/Lab/GiniClust/Rfunction/GiniClust_Clustering.R")
  Cluster.Results = GiniClust_Clustering(data.type = 'RNA-seq',ExprM.RawCounts.filter=data,GeneList.final=GeneList.final,eps=eps,MinPts=MinPts,out.folder=data_outdir,exprimentID=c)
  cell.cell.distance = Cluster.Results$cell_cell_dist
  c_membership = Cluster.Results$c_membership # group of each column(cells) 
  clustering_membership_r = Cluster.Results$clustering_membership_r # clusters of each cells with cells' name
  rare.cells.list.all = Cluster.Results$rare.cell
  
  #tSNE visualzation
  source("~/projects/def-cdesouza/Lab/GiniClust/Rfunction/GiniClust_tSNE.R")
  GiniClust_tSNE(data.type = 'RNA-seq',c_membership,cell.cell.distance,perplexity,out.folder=data_outdir,exprimentID=c)
  
  #check the clustering results
  table(c_membership)
  print(rare.cells.list.all)
  # final clustering of GiniClust
  write.csv(clustering_membership_r, file=paste0(data_outdir,"/clustersGiniCLust_",c,"_",args$name_dataset,".csv"))
  
}

print("DONE")

#======================================================================================================