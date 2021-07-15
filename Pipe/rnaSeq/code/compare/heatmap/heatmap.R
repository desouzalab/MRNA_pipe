#======================================================================================================
# COMPARE - heatmap
#======================================================================================================

#======================
# libraries
#======================
# For heatmap
suppressMessages(library(argparse))
suppressMessages(library(openxlsx))
suppressMessages(library(pheatmap))
suppressMessages(library(cowplot))
suppressMessages(library(data.table))


RSCRIPT <- "Rscript"

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

all_backspin_cluster_ssRNASeq_files <- list.files(args$backspin_cluster_input_directory, pattern = "*.csv*")
print(all_backspin_cluster_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files) & length(all_giniclust_cluster_ssRNASeq_files)==length(all_backspin_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){
    print(c)
    
    ### Create data frame
    # Read .csv file containing preprocessed data
    data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[c]),row.names=1)
    print("  ...read")
    data=na.omit(data)
    data=as.matrix(data)
    
    ### Load data
    # Read .xlsx file containing true cluster data
    if(args$name_dataset == "LaManno"){
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))
    apply_paste<-function(x){
      paste("X",x,sep="")
    } 
    true$GSM.ID<- sapply(true$GSM.ID, apply_paste)
    true=as.data.frame(setDT(true)[true$GSM.ID %chin% colnames(data)])[,2]
    
    }
    else{
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))
    true=as.data.frame(setDT(true)[true$GSM.ID %chin% colnames(data)])[,3]
    }
    # Read .csv file containing sc3 cluster data
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]))[,2]
    sc3=as.data.frame(sc3)
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    seurat=as.data.frame(seurat)
    # Read .csv file containing giniclust cluster data
    giniclust=read.csv(file.path(args$giniclust_cluster_input_directory, all_giniclust_cluster_ssRNASeq_files[c]))[,3]
    giniclust=as.data.frame(giniclust)
    # Read .csv file containing backspin cluster data
    backspin=read.csv(file.path(args$backspin_cluster_input_directory, all_backspin_cluster_ssRNASeq_files[c]))[,2]
    backspin=as.data.frame(backspin)
    #Extract 100 most variable genes
    #synthesise example data matrix
    if(args$name_dataset == "GSE74672"){
      data=data[rowSums(data)<6000 & rowSums(data)>200,]
      # NROW BEFORE 18548
    }


    print(nrow(data))
    mostVar <- function(data, n, i_want_most_var = TRUE) {
      data.var <- apply(data, 1, stats::var)
      data[order(data.var, decreasing = i_want_most_var)[1:n],] 
    }
    mv_100 <- mostVar(data, n = 100)
    ### Plots
    print(rowSums(mv_100))
    forhm = data.frame(TrueClusters=true,Seurat=seurat,SC3=sc3,GiniClust=giniclust,Backspin=backspin)

    labelsforhm=forhm[with(forhm,order(TrueClusters)),] # sort by true group
    rownames(labelsforhm)=colnames(mv_100)
    labelsmv_100 = mv_100[,rownames(true)] #sort in column (true group of cells)
    #HMTrueFirst=pheatmap(t(mv_100),cluster_rows=F,cluster_cols=F,scale="none",show_rownames=FALSE, show_colnames=FALSE, annotation_row=labelsforhm,fontsize=5,height=10,width=10) 
    data_plus=data[row.names(data) %in% row.names(mv_100),]
    cat(rowSums(data_plus),"this is plus :)")
    HMTrueFirst=pheatmap(t(data_plus),cluster_rows=F,cluster_cols=F,scale="row",show_rownames=FALSE, show_colnames=FALSE, annotation_row=labelsforhm,fontsize =5) 

    save_plot(paste0(outdir,"/heatmapTrueFirst_",c,"_",args$name_dataset,".png"),HMTrueFirst)
    print(rownames(data)[1:10])
    data_minus=data[!(row.names(data) %in% row.names(mv_100)),]

    #HMTrueFirst=pheatmap(t(mv_100),cluster_rows=F,cluster_cols=F,scale="none",show_rownames=FALSE, show_colnames=FALSE, annotation_row=labelsforhm,fontsize=5,height=10,width=10) 
    HMTrueFirst=pheatmap(t(data_minus),cluster_rows=F,cluster_cols=F,scale="none",show_rownames=FALSE, show_colnames=FALSE, annotation_row=labelsforhm,fontsize =6) 

    save_plot(paste0(outdir,"/heatmapTrueFirst_no_var_genes",c,"_",args$name_dataset,".png"),HMTrueFirst)
    
    print("  ...plot heatmapTrueFirst")


    rm(dat)
    rm(data)
    rm(data1)
    rm(true)
    rm(sc3)
    rm(seurat)
    rm(labelsforhm)
    rm(HMTrueFirst)
    rm(HMSeuratFirst)
  }
}
print("DONE")
#=====================================================================================================

​