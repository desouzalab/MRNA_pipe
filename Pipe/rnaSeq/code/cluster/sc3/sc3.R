#======================================================================================================
# CLUSTER - SC3
#======================================================================================================

#======================
# argparse library
#======================
suppressMessages(library(argparse))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(SC3))
suppressMessages(library(scater))
suppressMessages(library(cowplot))


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
  data=read.csv(file.path(args$input_directory, all_preprocessed_ssRNASeq_files[c]))
  print("  ...read")
  ### Set row names for the data frame. Exclude the first column from the data frame.
  row.names(data)=data[,1]

  ### Exclude the first column from the data frame.
  data=data[,-1]

  ### Cluster
  sce <- SingleCellExperiment(
      assays  = list(
          counts = as.matrix(data),
          normcounts = t(t(as.matrix(data))/colSums(as.matrix(data)))*1000000,
          logcounts = log2(as.matrix(data) + 1)))


  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

  sce <- runPCA(sce)

  ### Plot
  plot1=plotPCA(sce)
  save_plot(paste0(data_outdir,"/PCA_",c,"_",args$name_dataset,".pdf"),plot1)
  dev.off()
  print("  ...plot PCA")

  sce <- sc3_estimate_k(sce)
    optimal_K <- metadata(sce)$sc3$k_estimation

  sce <- sc3(sce, ks=optimal_K, biology=F)

  p_Data <- colData(sce)
    col_name <- paste("sc3_", optimal_K, "_clusters", sep='')
    sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]

  write.table(sc3OUTPUT, file=paste0(data_outdir,"/clustersSC3_",c,"_",args$name_dataset,".csv"), sep="\t", row.names=T)
  print("  ...export to .csv")

  rm(data)
  rm(plot1)
  rm(sce)
  rm(optimal_K)
  rm(p_data)
  rm(col_name)
  rm(sc3OUTPUT)

}

print("DONE")

#======================================================================================================
