#!/bin/bash

CLUSTERMETHOD=backspin
DATASET=zeisel
module load r/4.0.0

#sc3
#Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R 
--data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} 
--console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} 
--input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}
 --name_dataset ${DATASET}

#seurat
#Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R --data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/raw/${DATASET} --name_dataset ${DATASET}

#GINI
#Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R 
--data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} 
--console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} 
--input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}
 --name_dataset ${DATASET}
 --input_raw /home/emiliano/projects/def-cdesouza/Lab/data/raw/${DATASET}

#call spin
#sbatch --job-name="SPIN ${DATASET}" --output=/home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/cluster/backspin/${DATASET}/CALL_SPIN.out --export=DATASET=${DATASET} /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/cluster/call_spin.sh

#backspin cef cluster convert
Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R --data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --input_directory /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/cluster/backspin/${DATASET}  --name_dataset ${DATASET}

