
#======================================================================================================
# COMPARE - AIR+VM+Purity
#======================================================================================================

#======================
# libraries
#======================
# For AIR+VM+Purity
suppressMessages(library(argparse))
suppressMessages(library(openxlsx))
suppressMessages(library(mclust))
suppressMessages(library(clevr))
suppressMessages(library(NMF))
suppressMessages(library(ggplot2))
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

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.xlsx*")
print(all_true_cluster_ssRNASeq_files)

all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, pattern = "*.csv*")
print(all_sc3_cluster_ssRNASeq_files)

all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, pattern = "*.csv*")
print(all_seurat_cluster_ssRNASeq_files)

if (length(all_true_cluster_ssRNASeq_files)==length(all_sc3_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files)) {
  for (c in 1:length(all_true_cluster_ssRNASeq_files)){
    print(c)

    ### Load data
    true=read.xlsx(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]), sep="\t")
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]), sep="\t")
    print("  ...read")

    sc3=as.matrix(sc3)
    seurat=as.matrix(seurat)
    
    ### AIR
    AIRsc3=adjustedRandIndex(true,sc3)
    outputString <- "AIR SC3:"
    print(outputString)
    print(AIRsc3)

    AIRseurat=adjustedRandIndex(true,seurat)
    outputString <- "AIR SEURAT:"
    print(outputString)
    print(AIRseurat)

    AIRSeuratsc3=adjustedRandIndex(seurat,sc3)
    outputString <- "AIR SEURAT + SC3:"
    print(outputString)
    print(AIRSeuratsc3)


    ### V_measure
    Vsc3=v_measure(as.numeric(true),as.numeric(sc3))
    outputString <- "V MEASURE SC3:"
    print(outputString)
    print(Vsc3)

    Vsc3=v_measure(as.numeric(true),as.numeric(seurat))
    outputString <- "V MEASURE SEURAT:"
    print(outputString)
    print(Vsc3)


    ### purity (not needed)
    homSC3=purity(as.factor(true),as.factor(sc3))
    outputString <- "PURITY SC3:"
    print(outputString)
    print(homSC3)

    homSeu=purity(as.factor(true),as.factor(seurat))
    outputString <- "\nPURITY SEURAT:"
    print(outputString)
    print(homSeu)


    ### plot
    Method=c("Seurat","SC3")
    hom=c(0.4199931,0.3648039)
    clusters=c(19,33)

    hommat=data.frame(Method,hom,clusters)

    homplot=ggplot(hommat, aes(y=hom,x=clusters,color=Method))+geom_point(aes(color=Method))+theme(legend.position = "none")+expand_limits(x=c(0,50), y=c(0, 1))+ labs(x = "Number of Clusters", y = "Purity")+geom_text(aes(label=Method),hjust=0, vjust=2)+geom_vline(xintercept=7, linetype="dashed", color = "green")
    save_plot(paste0(outdir,"/purity_",c,"_",args$name_dataset,".pdf"),homplot)
    dev.off()
    print("  ...plot purity")


    rm(true)
    rm(sc3)
    rm(seurat)
    rm(AIRsc3)
    rm(Vsc3)
    rm(homSC3)
    rm(homSeu)
    rm(outputString)
    rm(Method)
    rm(hom)
    rm(clusters)
    rm(hommat)
    rm(homplot)
  }
}

print("DONE")
#======================================================================================================



