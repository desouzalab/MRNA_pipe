#!/bin/bash

module load r/4.0.0
VISUALIZEMETHOD=tSNE+PCA
DATASET=GSE98816+58
CLUSTERMETHOD=sc3
Rscript ../../rnaSeq/code/visualize/${VISUALIZEMETHOD}/${VISUALIZEMETHOD}.R --output_directory ../../rnaSeq/output/visualize/${VISUALIZEMETHOD}/${CLUSTERMETHOD}/${DATASET}/ --preprocessed_input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}/ --true_cluster_input_directory ~/projects/def-cdesouza/Lab/data/true/${DATASET}/ --cluster_input_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET}/ --name_dataset ${DATASET}


