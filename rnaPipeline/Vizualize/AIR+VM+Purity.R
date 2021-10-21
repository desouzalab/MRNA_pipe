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


#===========================================ANALYSIS===========================================#

print(c)
sink(paste0(outdir,"/",args$name_dataset,".txt"))
### Create data frame
# Read .csv file containing preprocessed data
data=read.csv(file.path(args$preprocessed_input_directory, all_preprocessed_ssRNASeq_files[1]),row.names=1)
print("  ...read")
data=na.omit(data)
data=as.matrix(data)

### Load data
# Read .xlsx file containing true cluster data
if(args$name_dataset == "LaManno" | args$name_dataset == "zeisel" ){
  true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[1]))
  apply_paste<-function(x){
    paste("X",x,sep="")
  } 
    
  true$GSM.ID<- sapply(true$GSM.ID, apply_paste)
  true=true[true$GSM.ID %in% colnames(data),]
  true=true[,2]
} else {
  true=read.csv(file.path(args$true_cluster_input_directory, all_true_cluster_ssRNASeq_files[1]))
  true=as.data.frame(setDT(true)[true$GSM.ID %chin% colnames(data)])[,2]
}
# Read .csv file containing sc3 cluster data

sc3=read.csv(file.path(args$sc3_cluster_input_directory, all_sc3_cluster_ssRNASeq_files[1]))[,2]
# Read .csv file containing seurat cluster data
seurat=read.csv(file.path(args$seurat_cluster_input_directory, all_seurat_cluster_ssRNASeq_files[1]))[,2]
giniclust=read.csv(file.path(args$giniclust_cluster_input_directory, all_giniclust_cluster_ssRNASeq_files[1]))[,3]
giniclust = as.numeric(factor(as.data.frame(giniclust)[,1])) 

backspin=read.csv(file.path(args$backspin_cluster_input_directory, all_backspin_cluster_ssRNASeq_files[1]))
print(head(backspin))
backspin=as.data.frame(backspin)
backspin=as.data.frame(backspin[backspin[,1] %in% colnames(data),])
#backspin=backspin[match(backspin[,1], colnames(data)),]

sc3=as.matrix(sc3)
seurat=as.matrix(seurat)
giniclust=as.matrix(giniclust)
backspin=as.matrix(backspin[,2])
print(length(backspin))
print(length(true))
### AIR
AIRsc3=adjustedRandIndex(true,sc3)
outputString <- "AIR SC3:"
print(outputString)
print(AIRsc3)

AIRseurat=adjustedRandIndex(true,seurat)
outputString <- "AIR SEURAT:"
print(outputString)
print(AIRseurat)

AIRginiclust=adjustedRandIndex(true,giniclust)
outputString <- "AIR giniclust:"
print(outputString)
print(AIRginiclust)



AIRbackspin=adjustedRandIndex(true,backspin)
outputString <- "AIR backspin:"
print(outputString)
print(AIRbackspin)


#V_measure
vSC3=v_measure(as.numeric(as.factor(true)),as.numeric(sc3))
outputString <- "v_measure SC3:"
print(outputString)
print(vSC3)

vSeu=v_measure(as.numeric(as.factor(true)),as.numeric(seurat))
outputString <- "\nv_measure SEURAT:"
print(outputString)
print(vSeu)

vSpin=v_measure(as.numeric(as.factor(true)),as.numeric(backspin))
outputString <- "\nv_measure BACKSPIN:"
print(outputString)
print(vSpin)

vGIniclust=v_measure(as.numeric(as.factor(true)),as.numeric(giniclust))
outputString <- "\nv_measure GINICLUST:"
print(outputString)
print(vSeu)


### purity (not needed)
homSC3=purity(as.factor(true),as.factor(sc3))
outputString <- "PURITY SC3:"
print(outputString)
print(homSC3)

homSeu=purity(as.factor(true),as.factor(seurat))
outputString <- "\nPURITY SEURAT:"
print(outputString)
print(homSeu)

homSpin=purity(as.factor(true),as.factor(backspin))
outputString <- "\nPURITY backspin:"
print(outputString)
print(homSpin)

homGIniclust=purity(as.factor(true),as.factor(giniclust))
outputString <- "\nPURITY giniclust:"
print(outputString)
print(homGIniclust)



### plot purity
Method=c("Seurat","SC3","Backspin","Giniclust")
hom=c(homSeu,homSC3,homSpin,homGIniclust)
clusters=c(length(levels(as.factor(seurat))),length(levels(as.factor(sc3))),length(levels(as.factor(backspin))),length(levels(as.factor(giniclust))))

hommat=data.frame(Method,hom,clusters)

homplot=ggplot(hommat, aes(y=hom,x=clusters,color=Method))+geom_point(aes(color=Method))+theme(legend.position = "none")+expand_limits(x=c(0,50), y=c(0, 1))+ labs(x = "Number of Clusters", y = "Purity")+geom_text(aes(label=Method),hjust=0, vjust=2)+geom_vline(xintercept=length(levels(as.factor(true))), linetype="dashed", color = "green")
save_plot(paste0(outdir,"/purity_",args$name_dataset,".pdf"),homplot)
dev.off()
print("  ...plot purity")

### plot V-measure
Method=c("Seurat","SC3","Backspin","Giniclust")
hom=c(vSeu,vSC3,vSpin,vGIniclust)
clusters=c(length(levels(as.factor(seurat))),length(levels(as.factor(sc3))),length(levels(as.factor(backspin))),length(levels(as.factor(giniclust))))

hommat=data.frame(Method,hom,clusters)

homplot=ggplot(hommat, aes(y=hom,x=clusters,color=Method))+geom_point(aes(color=Method))+theme(legend.position = "none")+expand_limits(x=c(0,50), y=c(0, 1))+ labs(x = "Number of Clusters", y = "V index")+geom_text(aes(label=Method),hjust=0, vjust=2)+geom_vline(xintercept=length(levels(as.factor(true))), linetype="dashed", color = "green")
save_plot(paste0(outdir,"/v_measure_",args$name_dataset,".pdf"),homplot)
dev.off()
print("  ...plot v measure")
sink()
rm(list = ls())

print("DONE")