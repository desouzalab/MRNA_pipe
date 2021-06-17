#======================================================================================================
# CLUSTER - SEURAT
#======================================================================================================

#======================
# libraries
#======================
suppressMessages(library(argparse))
library(Seurat)
library(cowplot)

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
  print(c)
  ### Create data frame
  dat=read.csv(file.path(args$input_directory, all_preprocessed_ssRNASeq_files[c]))
  data=dat
  print("  ...read")
  print(head(data))
  ### Set row names for the data frame. Exclude the first column from the data frame.
  row.names(data)=data[,1]

  ### Exclude the first column from the data frame.
  data=data[,-1]

  ### Run Seurat
  pbmc <- CreateSeuratObject(counts=data, project="pbmc3k", min.cells=3, min.features=200)

  outputString <- "PBMC:"
  print(outputString)
  pbmc

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


  # Determine the ‘dimensionality’ of the dataset
  plot3=ElbowPlot(object=pbmc, ndims=50)
  save_plot(paste0(data_outdir,"/elbow_",c,"_",args$name_dataset,".pdf"),plot3)
  dev.off()
  print("  ...plot elbow")


  # Cluster the cells
  pbmc <- FindNeighbors(pbmc, dims = 1:7)

  # Not being used to find clusters
  # CH=matrix(c(paste("CH-",seq(0.06,0.19,by=0.01),sep=""),rep(0,14)),ncol=2)
  # kclust0=0
  # for (i in 6:14){
  # pbmc1=pbmc
  # pbmc1 <- FindClusters(pbmc1, resolution = i/100)
  # clu=as.matrix(Idents(pbmc1))
  # clu=as.integer(clu[,1])
  # CH[i,2]=calinhara(t(data),clustering=clu,cn=max(clu))
  # }

  # if(CH>CH0){
  # CH0<-CH
  # kclust=
  # res=i/100} 
  # Not being used up to this point


  # Clusters
  pbmc <- FindClusters(pbmc, resolution = 1)
  write.table(Idents(pbmc), file=paste0(data_outdir,"/clustersSeurat_",c,"_",args$name_dataset,".csv"),sep="\t", row.names=T)
  print("  ...export to .csv")

  # Plots
  plot5=DimPlot(pbmc, reduction = "pca")
  save_plot(paste0(data_outdir,"/PCA_",c,"_",args$name_dataset,".pdf"),plot5)
  print("  ...plot PCA")

  pbmc <- RunUMAP(pbmc, dims = 1:10)
  plot6=DimPlot(pbmc, reduction = "umap")
  save_plot(paste0(data_outdir,"/UMAP_",c,"_",args$name_dataset,".pdf"),plot6)
  print("  ...plot UMAP")

  
  pbmc=RunTSNE(pbmc)
  plot7=DimPlot(pbmc, reduction = "tsne")
  save_plot(paste0(data_outdir,"/TSNE_",c,"_",args$name_dataset,".pdf"),plot7)
  print("  ...plot TSNE")
  dev.off()

  rm(dat)
  rm(data)
  rm(outputString)
  rm(pbmc)
  rm(top20)
  rm(plot1)
  rm(plot2)
  rm(combpo)
  rm(plot3)
  rm(plot5)
  rm(plot6)
  rm(plot7)
}

print("done")

#======================================================================================================
