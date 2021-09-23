
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
library(data.table)


RSCRIPT <- "Rscript"

#======================
# arguments
#======================

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--preprocessed_input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--true_cluster_input_directory", default="None", type="character", help="path to directory containing true cluster ssRNASeq data")
parser$add_argument("--sc3_cluster_input_directory", default="None", type="character", help="path to directory containing sc3 clustered ssRNASeq data")
parser$add_argument("--seurat_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")
parser$add_argument("--backspin_cluster_input_directory", default="None", type="character", help="path to directory containing backspin clustered ssRNASeq data")
parser$add_argument("--giniclust_cluster_input_directory", default="None", type="character", help="path to directory containing giniclust clustered ssRNASeq data")

parser$add_argument("--output_directory", type="character", help="Path to the output directory")

args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv*")
print(all_true_cluster_ssRNASeq_files)

all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, pattern = "*.csv*")
print(all_sc3_cluster_ssRNASeq_files)

all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, pattern = "*.csv*")
print(all_seurat_cluster_ssRNASeq_files)

all_giniclust_cluster_ssRNASeq_files <- list.files(args$giniclust_cluster_input_directory, pattern = "*.csv*")
print(all_giniclust_cluster_ssRNASeq_files)

all_backspin_cluster_ssRNASeq_files <- list.files(args$backspin_cluster_input_directory, pattern = "*.csv")
print(all_backspin_cluster_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files) & length(all_giniclust_cluster_ssRNASeq_files)==length(all_backspin_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){

    sink(paste0(outdir,"/",args$name_dataset,".txt"))
    ### Create data frame
    # Read .csv file containing preprocessed data
    data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[c]),row.names=1)
    print("  ...read")
    data=na.omit(data)
    data=as.matrix(data)
    ### Load data
    # Read .xlsx file containing true cluster data
    if(args$name_dataset == "LaManno" | args$name_dataset == "zeisel" ){
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))
    apply_paste<-function(x){
      paste("X",x,sep="")
    } 
      
    true$GSM.ID<- sapply(true$GSM.ID, apply_paste)
    true=true[true$GSM.ID %in% colnames(data),]
    true=true[,2]
    true=as.factor(true)
    }
    else{
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))
    true=as.data.frame(setDT(true)[true$GSM.ID %chin% colnames(data)])[,3]
    true=as.factor(true)
    }
    # Read .csv file containing sc3 cluster data
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]))[,2]
    #sc3=as.factor(sc3)
    
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    #seurat=as.factor(seurat)
    # Read .csv file containing giniclust cluster data
    giniclust=read.csv(file.path(args$giniclust_cluster_input_directory, all_giniclust_cluster_ssRNASeq_files[c]))[,3]
    giniclust=as.factor(giniclust)
    # Read .csv file containing backspin cluster data
    backspin=read.csv(file.path(args$backspin_cluster_input_directory, all_backspin_cluster_ssRNASeq_files[c]))[,2]
    #backspin=as.factor(backspin)
    #Extract 100 most variable genes
    #synthesise example data matrix

    print("  ...read")
    data = na.omit(data)
    data=as.matrix(data)
    #=====================FORMAT DATA=====================#
    set.seed(123)
    print(head(sc3))
    print((length(levels(true))))
    
    tsnepca=Rtsne(X=t(data), dims=2, perplexity=30, theta=0, check_duplicates=F, pca=TRUE, partial_pca=FALSE, max_iter=1000, verbose=T, is_distance=FALSE, Y_init=NULL, pca_center=TRUE, pca_scale=F, normalize=F) 
    tsneX=tsnepca$Y[,1]
    tsneY=tsnepca$Y[,2]

    hommat_sc3=data.frame(tsneX,tsneY,true,sc3)
    hommat_spin=data.frame(tsneX,tsneY,true,backspin)
    hommat_seurat=data.frame(tsneX,tsneY,true,seurat)
    hommat_ginilust=data.frame(tsneX,tsneY,true,giniclust)
    
    #=====================PLOT DATA=====================#
    tsnepca=ggplot(hommat_sc3, aes(y=tsneY,x=tsneX,color=sc3) ) + geom_point(aes(shape=true),size=1) + scale_shape_manual(values=seq(0,length(levels(true)))) + theme( legend.key.size=unit(.1,"inches"))  + scale_color_gradientn(colours = rainbow(nlevels(as.factor(sc3))))
    save_plot(paste0(outdir,"/TSNE+PCA_true_and_method+sc3","_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot SC3 ")
    dev.off()

    tsnepca=ggplot(hommat_spin, aes(y=tsneY,x=tsneX,color=backspin) ) + geom_point(aes(shape=true),size=1) + scale_shape_manual(values=seq(0,length(levels(true)))) + theme( legend.key.size=unit(.1,"inches")) + scale_color_gradientn(colours = rainbow(nlevels(as.factor(backspin))))
    save_plot(paste0(outdir,"/TSNE+PCA_true_and_method+spin","_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA black and white")
    dev.off()

    tsnepca=ggplot(hommat_ginilust, aes(y=tsneY,x=tsneX,color=giniclust) ) + geom_point(aes(shape=true),size=1) + scale_shape_manual(values=seq(0,length(levels(true)))) + theme( legend.key.size=unit(.1,"inches")) 
    save_plot(paste0(outdir,"/TSNE+PCA_true_and_method+giniclust","_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA black and white")
    dev.off()

    tsnepca=ggplot(hommat_seurat, aes(y=tsneY,x=tsneX,color=seurat) ) + geom_point(aes(shape=true),size=1) + scale_shape_manual(values=seq(0,length(levels(true)))) + theme( legend.key.size=unit(.1,"inches")) + scale_color_gradientn(colours = rainbow(nlevels(as.factor(seurat))))
    save_plot(paste0(outdir,"/TSNE+PCA_true_and_method+seurat","_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA black and white")
    dev.off()
    
    #=====================CLEAR MEMORY=====================#
    rm(list = ls())
    sink()
  }
}
print("DONE")