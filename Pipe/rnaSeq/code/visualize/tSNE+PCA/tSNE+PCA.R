
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
    ### Create data frame
    # Read .csv file containing preprocessed data
    data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[c]))
    print("  ...read")

    ### Set row names for the data frame. Exclude the first column from the data frame.
    row.names(data)=data[,1]

    ### Exclude the first column from the data frame.
    data=data[,-1]

    data=as.matrix(data)
    set.seed(123)

    tsnepca=Rtsne(X=t(data), dims=2, perplexity=30, theta=0, check_duplicates=F, pca=TRUE, partial_pca=FALSE, max_iter=1000, verbose=F, is_distance=FALSE, Y_init=NULL, pca_center=TRUE, pca_scale=F, normalize=F) 
    tsneX=tsnepca$Y[,1]
    tsneY=tsnepca$Y[,2]

    # colNames=F --> First row of data will not be used as column names. (If TRUE, the first row of data is used as column names)
    TrueClusters=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))[,3]
    
    # Select Row 2 and exclude Column 1 from the data frame.
    TrueClusters=as.factor(TrueClusters)

    clusters=read.csv(file.path(args$cluster_input_directory,all_clustered_ssRNASeq_files[c]))[,2]
    clusters=t(as.vector(clusters))

    hommat=data.frame(tsneX,tsneY,TrueClusters,clusters)
    Method=c("Seurat","SC3")
    #tsnepca=ggplot(hommat, aes(y=tsneY,x=tsneX))+theme(legend.position = "none")+expand_limits(x=c(0,50), y=c(0, 1))+ labs(x = "Number of Clusters", y = "Purity")+geom_text(aes(label=Method),hjust=0, vjust=2)+geom_vline(xintercept=7, linetype="dashed", color = "green")
    #save_plot(paste0(outdir,"/TSNE+PCA_",c,"_",args$name_dataset,".pdf"),tsnepca)
    #print("  ...plot tSNE+PCA")
    #dev.off()

    tsnepca = ggplot(hommat, aes(y=tsneY,x=tsneX)) + geom_point(aes(shape=as.factor(TrueClusters)),size=1) + scale_shape(solid = TRUE)  
    save_plot(paste0(outdir,"/TSNE+PCA_Black_",c,"_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA black and white")


    tsnepca = ggplot(hommat, aes(y=tsneY,x=tsneX, color=as.factor(TrueClusters))) + geom_point(aes(shape=as.factor(clusters)),size=1) + scale_shape(solid = TRUE)
    save_plot(paste0(outdir,"/TSNE+PCA_Colour_",c,"_",args$name_dataset,".pdf"),tsnepca)
    print("  ...plot tSNE+PCA colour")

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

