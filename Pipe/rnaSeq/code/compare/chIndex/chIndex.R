#======================================================================================================
# COMPARE - CH Index
#======================================================================================================

#======================
# libraries
#======================
# For CH Index
suppressMessages(library(argparse))
suppressMessages(library(openxlsx))
suppressMessages(library(fpc))

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


    CHtrue=calinhara(t(data),clustering=true,cn=max(true))
    outputString <- "CH TRUE:"
    print(outputString)
    print(CHtrue)

    CHSC3=calinhara(t(data),clustering=sc3,cn=max(sc3))
    outputString <- "CH SC3:"
    print(outputString)
    print(CHSC3)

    CHSeurat=calinhara(t(data),clustering=seurat,cn=max(seurat))
    outputString <- "CH SEURAT:"
    print(outputString)
    print(CHSeurat)


  
    rm(data)
    rm(true)
    rm(sc3)
    rm(seurat)
    rm(CHtrue)
    rm(CHSC3)
    rm(CHSeurat)
    rm(outputString)

  }
}

print("DONE")
#======================================================================================================



