
#======================================================================================================
# VISUALIZATION - tSNE
#======================================================================================================

#======================
# libraries
#======================
# For tSNE
suppressMessages(library(argparse))
suppressMessages(library(M3C))
suppressMessages(library(cowplot))

RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--preprocessed_input_directory", default="None", type="character", help="path to directory containing preprocessed ssRNASeq data")
parser$add_argument("--trueCluster_input_directory", default="None", type="character", help="path to directory containing true cluster ssRNASeq data")
parser$add_argument("--cluster_input_directory", default="None", type="character", help="path to directory containing clustered ssRNASeq data")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")

args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)

all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)

all_clustered_ssRNASeq_files <- list.files(args$cluster_input_directory, pattern = "*.csv*")
print(all_clustered_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_clustered_ssRNASeq_files)) {
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

    ### TSNE
    # Read .csv file containing clusters
    
    labels=read.csv(file.path(args$cluster_input_directory, all_clustered_ssRNASeq_files[c]))[,2]
    print(labels)
    labels=as.matrix(labels)
    # Plot tSNE
    tsnePlot=tsne(data,labels=as.factor(labels),perplex=30 ,seed=123,dotsize=0.5,axistextsize=12, legendtextsize=2)
    save_plot(paste0(outdir,"/tSNE_",c,"_",args$name_dataset,".pdf"),tsnePlot)
    dev.off()
    print("  ...tSNE")
  }
}

print("DONE")
#======================================================================================================
