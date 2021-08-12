#======================================================================================================
# CLUSTER - SEURAT
#======================================================================================================

#======================
# libraries
#======================
suppressMessages(library(argparse))
library(Seurat)
library(cowplot)
library(data.table)
RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--data_output_directory", type="character", help="Path to the cluster data and plot output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE, recursive=TRUE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)


for (c in 1:length(all_preprocessed_ssRNASeq_files)){
  cat(c)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_preprocessed_ssRNASeq_files[c]),row.names=1)
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

  # Identification of highly variable features (feature selection)
  pbmc <- FindVariableFeatures(pbmc, selection.method="vst", nfeatures=2000)

  # Identify the 20 most highly variable genes
  top20 <- head(VariableFeatures(pbmc), 20)
  
  outputString <- "20 MOST HIGHLY VARIABLE GENES:"
  print(outputString)
  print(top20)
  
  # plot variable features with and without labels (give user an option to do this)
  plot1 <- VariableFeaturePlot(pbmc)
  plot2 <- LabelPoints(plot=plot1, points=top20, repel=TRUE)
  combpo=CombinePlots(plots=list(plot1, plot2))
  
  # save plot (give user an option to do this)
  pdf(paste0(data_outdir,"/variableFeaturesWithoutLabels_",c,"_",args$name_dataset,".pdf"))
  plot1 
  dev.off()
  print("  ...plot variableFeaturesWithoutLabels")
  
  pdf(paste0(data_outdir,"/variableFeaturesWithLabels_",c,"_",args$name_dataset,".pdf"))
  plot2
  dev.off()
  print("  ...plot variableFeaturesWithLabels")
  
  
  # Scaling the data
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  
  # Perform linear dimensional reduction
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  
  # display PCA results (OPTIONAL)
  outputString <- "PCA RESULTS:"
  print(outputString)
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  
  
  # Determine the ¡®dimensionality¡¯ of the dataset
  plot3=ElbowPlot(object=pbmc, ndims=50)
  save_plot(paste0(data_outdir,"/elbow_",c,"_",args$name_dataset,".pdf"),plot3)
  dev.off()
  print("  ...plot elbow")
  
  
  # Cluster the cells
  pbmc <- FindNeighbors(pbmc, dims = 1:30)
  
  # Clusters
  pbmc <- FindClusters(pbmc, resolution = 1)
  print(pbmc)
  write.csv(Idents(pbmc), file=paste0(data_outdir,"/seurat_",args$name_dataset,".csv"))
  print("  ...export to .csv")
  
  # Plots
  plot5=DimPlot(pbmc, reduction = "pca")
  save_plot(paste0(data_outdir,"/PCA_",c,"_",args$name_dataset,".pdf"),plot5)
  print("  ...plot PCA")
  
  pbmc <- RunUMAP(pbmc, dims = 1:30)
  plot6=DimPlot(pbmc, reduction = "umap")
  save_plot(paste0(data_outdir,"/UMAP_",c,"_",args$name_dataset,".pdf"),plot6)
  print("  ...plot UMAP")

  pbmc=RunTSNE(pbmc)
  plot7=DimPlot(pbmc, reduction = "tsne")
  save_plot(paste0(data_outdir,"/TSNE_",c,"_",args$name_dataset,".pdf"),plot7)
  print("  ...plot TSNE")
  dev.off()
  
  rm(list = ls())

}

print("done")

#======================================================================================================