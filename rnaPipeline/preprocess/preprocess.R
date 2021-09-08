#======================================================================================================
# PREPROCESSNG
#======================================================================================================

#======================
# libraries
#======================
suppressMessages(library(argparse))
suppressMessages(library(argparse))
library(dplyr)
library(Seurat)
library(patchwork)
suppressMessages(library(cowplot))
RSCRIPT <- "Rscript"


#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing raw ssRNASeq data")
parser$add_argument("--data_output_directory", type="character", help="Path to the preprocessed data output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE)

all_raw_ssRNASeq_files <- list.files(args$input_directory, pattern="*.csv*")
print(all_raw_ssRNASeq_files)
print("hi")
for (c in 1:length(all_raw_ssRNASeq_files)){
  cat(c)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_raw_ssRNASeq_files[c]),row.names=1)
  print("  ...read")
  data=na.omit(data)
  #print(head(data[1:10]))
  pbmc <- CreateSeuratObject(counts = data, project = "data3k", min.cells = 3, min.features = 100)
  # focus on MT features: Low-quality / dying cells often exhibit extensive mitochondrial contamination
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  Idents(object=pbmc) <- "data3k"
  #Cutoff Plots
  vlnplot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  save_plot(paste0(args$console_output_directory,"/VlnPlot",c,"_",args$name_dataset,".pdf"),vlnplot)
  print("  ...plot tSNE+PCA colour")
  dev.off()

  # FeatureScatter to visualize feature-feature relationships
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  save_plot(paste0(args$console_output_directory,"/FeatureScatter",c,"_",args$name_dataset,".pdf"),plot1)
  print("  ...plot tSNE+PCA colour")
  dev.off()
  save_plot(paste0(args$console_output_directory,"/FeatureScatter_2",c,"_",args$name_dataset,".pdf"),plot2)
  print("  ...plot tSNE+PCA colour")
  dev.off()

  # Normalize the data
  pbmc <- NormalizeData(pbmc, normalization.method="LogNormalize", scale.factor=10000)
  outFilename <- paste0(data_outdir,"/",args$name_dataset,"_pre",".csv")
  write.csv(as.matrix(GetAssayData(pbmc, slot = "counts")), file=outFilename,row.names=TRUE)
  print("  ...export to .csv")

  rm(list = ls())

}

print("DONE")

#======================================================================================================

