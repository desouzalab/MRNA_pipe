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

for (c in 1:length(all_raw_ssRNASeq_files)){
  cat(c)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_raw_ssRNASeq_files[c]))
  print("  ...read")
  seurat_object <- CreateSeuratObject(counts = data, project = "data3k", min.cells = 3, min.features = 200)
  seurat_object
  # focus on MT features: Low-quality / dying cells often exhibit extensive mitochondrial contamination
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  vlnplot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  save_plot(paste0(args$console_output_directory,"/VlnPlot",c,"_",args$name_dataset,".pdf"),vlnplot)
  print("  ...plot tSNE+PCA colour")
  dev.off()
  # FeatureScatter to visualize feature-feature relationships
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  save_plot(paste0(args$console_output_directory,"/FeatureScatter",c,"_",args$name_dataset,".pdf"),plot1)
  print("  ...plot tSNE+PCA colour")
  dev.off()
  save_plot(paste0(args$console_output_directory,"/FeatureScatter_2",c,"_",args$name_dataset,".pdf"),plot2)
  print("  ...plot tSNE+PCA colour")
  dev.off()
  # based on figures, filtering (choose the threshold based on plots)
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  #Normalization
  pbmc <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)


  ### Exclude records where rowsum is less than 50
  data=data[!rowSums(data)<50,]

  ### Exclude records that have less than or equal to 864 zero's
  data=data[!apply(data==0, 1, sum) <= 864, ]
  print("  ...preprocess")

  ### Export preprocessed data frame to CSV file
  outFilename <- paste0(data_outdir,"/preprocessed_",c,"_",args$name_dataset,".csv")
  write.csv(data, file=outFilename,row.names=TRUE)
  print("  ...export to .csv")

  rm(dat)
  rm(data)
  rm(outFilename)
}

print("DONE")

#======================================================================================================

