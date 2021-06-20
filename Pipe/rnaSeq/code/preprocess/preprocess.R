#======================================================================================================
# PREPROCESSNG
#======================================================================================================

#======================
# libraries
#======================
suppressMessages(library(argparse))

RSCRIPT <- "Rscript"


#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing raw ssRNASeq data")
parser$add_argument("--data_output_directory", type="character", help="Path to the preprocessed data output directory")
parser$add_argument("--console_output_directory", type="character", help="Path to the .out file output directory")

args <- parser$parse_args()
print(args)

data_outdir <- args$data_output_directory
dir.create(file.path(data_outdir), showWarnings=FALSE)

console_outdir <- args$console_output_directory
dir.create(file.path(console_outdir), showWarnings=FALSE)

all_raw_ssRNASeq_files <- list.files(args$input_directory, pattern="*.csv*")
print(all_raw_ssRNASeq_files)

for (c in 1:length(all_raw_ssRNASeq_files)){
  cat(c)
  ### Create data frame
  data=read.csv(file.path(args$input_directory, all_raw_ssRNASeq_files[c]))
  print("  ...read")
  print(head(data))
  ### Set row names for the data frame. 
  row.names(data)=data[,1]

  ### Exclude the first column from the data frame.
  data=data[,-1]

  ### Exclude records where rowsum is less than 50
  data=data[!rowSums(data)<50,]

  ### Exclude records that have less than or equal to 864 zero's
  data=data[!apply(data==0, 1, sum) <= 864, ]
  print("  ...preprocess")

  ### Export preprocessed data frame to CSV file
  outFilename <- paste0(data_outdir,"/preprocessed_",c,"_",args$name_dataset,".csv")
  write.csv(data, file=outFilename,row.names=TRUE)
  print("  ...export to .csv")

  rm(dat)
  rm(data)
  rm(outFilename)
}

print("DONE")

#======================================================================================================

