#!/bin/bash

DATASET=GSE98816+58
module load r/4.0.0

Rscript /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/preprocess/Cef_convert.R --data_output_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --console_output_directory ../../rnaSeq/output/preprocess/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --name_dataset ${DATASET}
