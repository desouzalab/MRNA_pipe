#======================================================================================================
# COMPARE - Up and down regulated genes
#======================================================================================================

#======================
# libraries
#======================
# For up and down regulated genes
suppressMessages(library(argparse))
suppressMessages(library(openxlsx))
suppressMessages(library(tis))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(StatMeasures))

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

all_true_cluster_ssRNASeq_files <- list.files(args$true_cluster_input_directory, pattern = "*.csv*")
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
    true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[c]))[,3]
    # Read .csv file containing sc3 cluster data
    sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[c]))[,2]
    # Read .csv file containing seurat cluster data
    seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[c]))[,2]
    print("  ...read")
    sc3=as.matrix(sc3)
    seurat=as.matrix(seurat)

    ### Boxplot for True

    # add true to the last row of data
    data1=data
    data1[nrow(data1)+1,]<-true
    data1=t(data1)
    data1=data.frame(data1)

    # obtain mean of data based on the last row
    data1=cbind(as.matrix(t(data)),as.numeric(true))
    data1=data.frame(data1)

    avgtrue=aggregate(data1[,-13365], by=list(data1$V13365), FUN=mean)

    OutVals = boxplot(avgtrue[1,])$out
    which(avgtrue[1,] %in% OutVals)
    write.table(avgtrue, file=paste0(outdir,"/AverageTrue_",c,"_",args$name_dataset,".csv"),sep="\t", row.names=F)

    ### find first 5 top outliers in each cluster
    outputString <- "FIRST 5 TOP OUTLIERS IN EACH CLUSTER:"
    print(outputString)
    print(head(sort(avgtrue[1,],decreasing=F),n=5))
    print(head(sort(avgtrue[2,],decreasing=F),n=5))
    print(head(sort(avgtrue[3,],decreasing=F),n=5))
    print(head(sort(avgtrue[4,],decreasing=F),n=5))
    print(head(sort(avgtrue[5,],decreasing=F),n=5))
    print(head(sort(avgtrue[6,],decreasing=F),n=5))
    print(head(sort(avgtrue[7,],decreasing=F),n=5))

    ### Outliers
    avg1=avgtrue[1,]
    length(avg1[avg1>0.41,])
    out=outliers(avgtrue) 


    ### Boxplot for true

    avgtrue=data.frame(t(avgtrue))
    ClusTr=c(paste("C",seq(1,7,by=1),sep=""))

    pdf(file=paste0(outdir,"/BoxPlotTrue_",c,"_",args$name_dataset,".pdf"))
    boxplot(avgtrue, outcol="#FF33CC",names=ClusTr)
    dev.off()


    ### Boxplot for SC3

    #add true to the last row of data
    data1=data
    data1[nrow(data1)+1,]<-t(sc3)
    data1=t(data1)
    data1=data.frame(data1)

    # obtain mean of data based on the last row
    data1=cbind(as.matrix(t(data)),as.numeric(t(sc3)))
    data1=data.frame(data1)
    avgsc3=aggregate(data1[,-13365], by=list(data1$V13365), FUN=mean)
    avgsc3=data.frame(t(avgsc3))[-1,]

    ClusSc=c(paste("C",seq(1,33,by=1),sep=""))
    pdf(file=paste0(outdir,"/BoxPlotSC3_",c,"_",args$name_dataset,".pdf"))
    boxplot(avgsc3, outcol="green",names=ClusSc)
    dev.off()


    ### Boxplot for seurat

    # add true to the last row of data
    data1=data
    data1[nrow(data1)+1,]<-t(seurat)
    data1=t(data1)
    data1=data.frame(data1)

    #obtain mean of data based on the last row
    data1=cbind(as.matrix(t(data)),as.numeric(t(seurat)))
    data1=data.frame(data1)
    avgseu=aggregate(data1[,-13365], by=list(data1$V13365), FUN=mean)
    avgseu=data.frame(t(avgseu))[-1,]

    ClusSe=c(paste("C",seq(1,19,by=1),sep=""))
    pdf(file=paste0(outdir,"/BoxPlotSeurat_",c,"_",args$name_dataset,".pdf"))
    boxplot(avgseu, outcol="blue",names=ClusSe)
    dev.off()

    rm(dat)
    rm(data)
    rm(data1)
    rm(true)
    rm(sc3)
    rm(seurat)
    rm(avgtrue)
    rm(OutVals)
    rm(outputString)
    rm(avg1)
    rm(out)
    rm(ClusTr)
    rm(avgsc3)
    rm(ClusSc)
    rm(avgseu)

  }
}

print("DONE")
#======================================================================================================



