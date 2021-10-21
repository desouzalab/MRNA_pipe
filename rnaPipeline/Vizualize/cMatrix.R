#======================================================================================================
# COMPARE - Confusion Matrix
#======================================================================================================

#======================
# libraries
#======================
# For confusion matrix
suppressMessages(library(argparse))
suppressMessages(library(openxlsx))
suppressMessages(library(descr))
suppressMessages(library(pheatmap))
suppressMessages(library(cowplot))
suppressMessages(library(data.table))

RSCRIPT <- "Rscript"

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--preprocessed_input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--true_cluster_input_directory", default="None", type="character", help="path to directory containing true cluster ssRNASeq data")
parser$add_argument("--sc3_cluster_input_directory", default="None", type="character", help="path to directory containing sc3 clustered ssRNASeq data")
parser$add_argument("--seurat_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")
parser$add_argument("--backspin_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")
parser$add_argument("--giniclust_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")

parser$add_argument("--output_directory", type="character", help="Path to the output directory")

args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern =  paste0(args$name_dataset,"_pre.csv*"))
print(all_preprocessed_ssRNASeq_files)

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv")
print(all_true_cluster_ssRNASeq_files)

all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, paste0("*sc3_",args$name_dataset,".csv*"))
print(all_sc3_cluster_ssRNASeq_files)

all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, paste0("*seurat_",args$name_dataset,".csv*"))
print(all_seurat_cluster_ssRNASeq_files)

all_giniclust_cluster_ssRNASeq_files <- list.files(args$giniclust_cluster_input_directory, pattern = paste0("*giniclust_",args$name_dataset,".csv*"))
print(all_giniclust_cluster_ssRNASeq_files)

all_backspin_cluster_ssRNASeq_files <- list.files(args$backspin_cluster_input_directory, paste0("*backspin_",args$name_dataset,".csv*"))
print(all_backspin_cluster_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_true_cluster_ssRNASeq_files)==length(all_sc3_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){
    sink(paste0(outdir,"/",args$name_dataset,".txt"))
    ### Create data frame
    # Read .csv file containing preprocessed data
    data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[1]),row.names=1)
    print("  ...read")
    data=na.omit(data)
    data=as.matrix(data)
    if(args$name_dataset == "LaManno" | args$name_dataset == "zeisel" ){
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[1]))
    apply_paste<-function(x){
      paste("X",x,sep="")
    } 
      
    true$GSM.ID<- sapply(true$GSM.ID, apply_paste)
    true=true[true$GSM.ID %in% colnames(data),]
    true=true[,2]
    }
    else{
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[1]))
    true=as.data.frame(setDT(true)[true$GSM.ID %chin% colnames(data)])[,3]
    }
    #  Read .csv file containing sc3 cluster data
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[1]))[,2]
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[1]))[,2]
    giniclust=read.csv(file.path(args$giniclust_cluster_input_directory, all_giniclust_cluster_ssRNASeq_files[1]))[,3]
    backspin=read.csv(file.path(args$backspin_cluster_input_directory, all_backspin_cluster_ssRNASeq_files[1]))[,2]
    my_palette <- colorRampPalette(c("white","green"))(n = 299)
    
    Hmcstr=pheatmap( crosstab(true ,sc3,prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=levels(true), show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/Matrix_sc3_true_",args$name_dataset,".jpg"),Hmcstr)

    Hmsetr=pheatmap( crosstab(true ,backspin,prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=levels(true),show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/Matrix_backspin_true_",args$name_dataset,".jpg"),Hmsetr)
    dev.off()

    Hmsetr=pheatmap( crosstab(true ,seurat,prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=levels(true),show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/Matrix_seurat_true_",args$name_dataset,".jpg"),Hmsetr)
    dev.off()
    
    Hmsetr=pheatmap( crosstab(true ,giniclust,prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=levels(true),show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/Matri_Giniclust_true_",args$name_dataset,".jpg"),Hmsetr)
    dev.off()

    rm(list = ls())
    sink()
  }
}
print("DONE")