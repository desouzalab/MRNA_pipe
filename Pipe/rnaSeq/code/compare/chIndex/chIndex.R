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
parser$add_argument("--backspin_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")
parser$add_argument("--giniclust_cluster_input_directory", default="None", type="character", help="path to directory containing seurat clustered ssRNASeq data")

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
print("p0")
args <- parser$parse_args()
print(args)
print("pass0")
outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings=FALSE, recursive=TRUE)
print("pass1")
all_preprocessed_ssRNASeq_files <- list.files(args$preprocessed_input_directory, pattern = "*.csv*")
print(all_preprocessed_ssRNASeq_files)
print("pass2")
all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv*")
print(all_true_cluster_ssRNASeq_files)
print("pass3")
all_sc3_cluster_ssRNASeq_files <- list.files(args$sc3_cluster_input_directory, pattern = "*.csv*")
print(all_sc3_cluster_ssRNASeq_files)
print("pass4")
all_seurat_cluster_ssRNASeq_files <- list.files(args$seurat_cluster_input_directory, pattern = "*.csv*")
print(all_seurat_cluster_ssRNASeq_files)
print("pass5")
all_giniclust_cluster_ssRNASeq_files <- list.files(args$giniclust_cluster_input_directory, pattern = "*.csv*")
print(all_giniclust_cluster_ssRNASeq_files)
print("pass6")
all_backspin_cluster_ssRNASeq_files <- list.files(args$backspin_cluster_input_directory, pattern = "*.csv*")
print(all_backspin_cluster_ssRNASeq_files)

if (length(all_preprocessed_ssRNASeq_files)==length(all_true_cluster_ssRNASeq_files) & length(all_true_cluster_ssRNASeq_files)==length(all_sc3_cluster_ssRNASeq_files) & length(all_sc3_cluster_ssRNASeq_files)==length(all_seurat_cluster_ssRNASeq_files)) {
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
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    giniclust=read.csv(file.path(args$giniclust_cluster_input_directory, all_giniclust_cluster_ssRNASeq_files[c]))[,3]
    backspin=read.csv(file.path(args$backspin_cluster_input_directory, all_backspin_cluster_ssRNASeq_files[c]))[,2]
    print("  ...read")
    
    # Transform cluster results in proper structure
    cellnames = colnames(data) #cell names from preprocessed matrix
    print("pass1")
    true = as.data.frame(true)
    print("pass2")
    row.names(true)=cellnames
    print("pass3")
    sc3=as.data.frame(sc3)
    print("pass4")
    row.names(sc3)=cellnames
    print("pass5")
    seurat=as.data.frame(seurat)
    print("pass6")
    row.names(seurat)=cellnames
    print("pass7")
    giniclust=as.data.frame(giniclust)
    print("pass8")
    giniclust = as.numeric(factor(giniclust[,1])) #numeric it
    print("pass9")
    giniclust=as.data.frame(giniclust)
    print("pass10")
    row.names(giniclust)=cellnames 
    print("pass11")
    backspin=as.data.frame(backspin)
    print("pass12")
    row.names(backspin)=cellnames
    print("pass13")
    print(true)
    CHtrue=calinhara(t(data),clustering=as.numeric(as.factor(true)),cn=length(levels(as.factor(true))))
    outputString <- "CH TRUE:"
    print(outputString)
    print(CHtrue)
    
    CHSC3=calinhara(t(data),clustering=as.numeric(as.factor(sc3)),cn=max(sc3))
    outputString <- "CH SC3:"
    print(outputString)
    print(CHSC3)
    
    CHSeurat=calinhara(t(data),clustering=as.numeric(as.factor(seurat)),cn=max(seurat))
    outputString <- "CH SEURAT:"
    print(outputString)
    print(CHSeurat)
    
    CHGiniClust=calinhara(t(data),clustering=as.numeric(as.factor(giniclust)),cn=max(giniclust))
    outputString <- "CH GINICLUST:"
    print(outputString)
    print(CHGiniClust)
    
    CHBackSpin=calinhara(t(data),clustering=as.numeric(as.factor(backspin)),cn=max(BackSpin))
    outputString <- "CH BACKSPIN:"
    print(outputString)
    print(CHBackSpin)
    
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