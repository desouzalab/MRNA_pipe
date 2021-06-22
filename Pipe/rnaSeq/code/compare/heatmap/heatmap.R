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
parser$add_argument("--output_directory", type="character", help="Path to the output directory")

args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.xlsx*")
print(all_true_cluster_ssRNASeq_files)

all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, pattern = "*.csv*")
print(all_sc3_cluster_ssRNASeq_files)

all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, pattern = "*.csv*")
print(all_seurat_cluster_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_true_cluster_ssRNASeq_files)==length(all_sc3_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){
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
    data[data == 0] <- NA

    ### Load data
    # Read .xlsx file containing true cluster data
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))[,2]
    # Read .csv file containing sc3 cluster data
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]))[,2]
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    print("  ...read")

    sc3=as.matrix(sc3)
    seurat=as.matrix(seurat)

    ### Plots
    labelsforhm=data.frame(TrueClusters=as.factor(true),Seurat=as.factor(seurat),SC3=as.factor(sc3))
    rownames(labelsforhm) <- rownames(t(data))
    HMTrueFirst=pheatmap(t(data),cluster_rows=F,cluster_cols=F,scale="none",show_rownames=FALSE, show_colnames=FALSE, annotation_row=labelsforhm) 
    save_plot(paste0(outdir,"/heatmapTrueFirst_",c,"_",args$name_dataset,".jpg"),HMTrueFirst)
    print("  ...plot heatmapTrueFirst")

    labelsforhm=data.frame(Seurat=as.factor(seurat),TrueClusters=as.factor(t(true)),SC3=as.factor(sc3))
    rownames(labelsforhm) <- rownames(t(data))
    HMSeuratFirst=pheatmap(t(data),cluster_rows=F,cluster_cols=F,scale="none",show_rownames=FALSE, show_colnames=FALSE, annotation_row = labelsforhm) 
    save_plot(paste0(outdir,"/heatmapSeuratFirst_",c,"_",args$name_dataset,".jpg"),HMSeuratFirst)
    print("  ...plot heatmapSeuratFirst")
    dev.off()

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
#======================================================================================================



