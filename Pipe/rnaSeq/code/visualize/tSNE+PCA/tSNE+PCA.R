
#======================================================================================================
# VISUALIZATION - tSNE+PCA
#======================================================================================================

#======================
# libraries
#======================
# For tSNE+PCA
suppressMessages(library(argparse))
suppressMessages(library(Rtsne))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(openxlsx))
suppressMessages(library(M3C))



RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--preprocessed_input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--true_cluster_input_directory", default="None", type="character", help="path to directory containing true cluster ssRNASeq data")
parser$add_argument("--cluster_input_directory", default="None", type="character", help="path to directory containing clustered ssRNASeq data")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")

args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv*")
print(all_true_cluster_ssRNASeq_files)

all_clustered_ssRNASeq_files <- list.files(args$cluster_input_directory, pattern = "*.csv*")
print(all_clustered_ssRNASeq_files)


if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_true_cluster_ssRNASeq_files)==length(all_clustered_ssRNASeq_files)) {
  for (c in 1:length(all_preprocessed_ssRNASeq_files)){
    print(c)
    #=====================READ DATA=====================#
    data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[c]),row.names = 1)
    TrueClusters=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))[,2]
    TrueClusters=as.factor(TrueClusters)
    clusters=read.csv(file.path(args$cluster_input_directory,all_clustered_ssRNASeq_files[c]))[,2]
    clusters=as.factor(clusters)
    print("  ...read")

    data=as.matrix(data)
    #=====================FORMAT DATA=====================#
    set.seed(123)
    tsnepca=Rtsne(X=t(data), dims=2, perplexity=30, theta=0, check_duplicates=F, pca=TRUE, partial_pca=FALSE, max_iter=1000, verbose=F, is_distance=FALSE, Y_init=NULL, pca_center=TRUE, pca_scale=F, normalize=F) 
    tsneX=tsnepca$Y[,1]
    tsneY=tsnepca$Y[,2]
    hommat=data.frame(tsneX,tsneY,TrueClusters,clusters)
    Method=c("Seurat","SC3")
    #=====================PLOT DATA=====================#
    tsnepca=ggplot(hommat, aes(y=tsneY,x=tsneX, color=clusters)) + geom_point(aes(shape=TrueClusters),size=1) + scale_shape_manual(values=seq(0,length(labels(TrueClusters))))
    save_plot(paste0(outdir,"/TSNE+PCA_Colour_",c,"_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA colour")
    dev.off()

    tsnepca=ggplot(hommat, aes(y=tsneY,x=tsneX)) + geom_point(aes(shape=TrueClusters),size=1) + scale_shape_manual(values=seq(0,length(levels(TrueClusters))))
    save_plot(paste0(outdir,"/TSNE+PCA_Black_",c,"_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA black and white")
    dev.off()

    tsnePlot=tsne(data,labels=TrueClusters,perplex=30 ,seed=123,dotsize=1,axistextsize=12, legendtextsize=10)
    save_plot(paste0(outdir,"/tSNE_",c,"_",args$name_dataset,".pdf"),tsnePlot)
    dev.off()
    print("  ...tSNE")
    #=====================CLEAR MEMORY=====================#
    rm(dat)
    rm(data)
    rm(tsneX)
    rm(tsneY)
    rm(TrueClusters)
    rm(clusters)
    rm(hommat)
    rm(tsnepca)
  }
}

print("DONE")

#======================================================================================================

