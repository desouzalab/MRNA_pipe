#!/bin/bash

CLUSTERMETHOD=sc3
DATASET=GSE74672
module load r/4.0.0

Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R --data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --name_dataset ${DATASET}


#Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R --data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --name_dataset ${DATASET} --input_raw /home/emiliano/projects/def-cdesouza/Lab/data/raw/${DATASET}

