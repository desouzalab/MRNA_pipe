#!/bin/sh
cd /home/emiliano/projects/def-cdesouza/Lab 
mkdir data/raw/$DATA_NAME
mkdir data/preprocessed/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/cluster/sc3/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/cluster/seurat/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/compare/AIR+VM+Purity/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/compare/chIndex/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/compare/cMatrix/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/compare/heatmap/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/compare/regGenes/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/preprocess/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE/seurat/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE/sc3/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE+PCA/seurat/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE+PCA/sc3/$DATA_NAME

touch data/preprocessed/$DATA_NAME/$DATA_NAME.csv