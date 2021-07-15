#!/bin/bash

DATASET=LaManno
module load r/4.0.0

#rest
#Rscript /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/preprocess/preprocess.R --data_output_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --console_output_directory ../../rnaSeq/output/preprocess/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/raw/${DATASET} --name_dataset ${DATASET}

#backspin cef conv
Rscript /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/preprocess/Cef_convert.R --data_output_directory ~/projects/def-cdesouza/Lab/${DATASET} --console_output_directory ../../rnaSeq/output/preprocess/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --name_dataset ${DATASET}

