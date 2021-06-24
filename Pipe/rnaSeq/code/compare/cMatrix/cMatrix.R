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

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv*")
print(all_true_cluster_ssRNASeq_files)

all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, pattern = "*.csv*")
print(all_sc3_cluster_ssRNASeq_files)

all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, pattern = "*.csv*")
print(all_seurat_cluster_ssRNASeq_files)

if (length(all_true_cluster_ssRNASeq_files)==length(all_sc3_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){
    print(c)

    ### Load data
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))[,3]
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]))[,2]
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    print("  ...read")

    sc3=as.matrix(sc3)
    seurat=as.matrix(seurat)

    ### Plots
    my_palette <- colorRampPalette(c("white","green"))(n = 299)

    Hmcstr=pheatmap( crosstab(true ,sc3[,1],prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=c("oligos", "astrocytes","ependymal","microglia","vsm","endothelial","neurons"), show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/heatmapscrtrue_",c,"_",args$name_dataset,".jpg"),Hmcstr)
    print("  ...plot heatmapscrtrue")

    Hmsetr=pheatmap( crosstab(true ,seurat[,1],prop.r=T)$prop.r ,cluster_rows=F,cluster_cols=F,show_rownames=T,labels_row=c("oligos", "astrocytes","ependymal","microglia","vsm","endothelial","neurons"),show_colnames=T,color =my_palette) 
    save_plot(paste0(outdir,"/heatmapseutrue_",c,"_",args$name_dataset,".jpg"),Hmsetr)
    print("  ...plot heatmapseutrue")
    dev.off()

    rm(true)
    rm(sc3)
    rm(seurat)
    rm(my_palette)
    rm(Hmcstr)
    rm(Hmsetr)
  }
}

print("DONE")
#======================================================================================================



