#!/bin/sh
cd /home/emiliano/projects/def-cdesouza/Lab 
DATA_NAME=GSE74672
New_method=giniclust
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

mkdir MRNA_pipe/Pipe/rnaSeq/output/cluster/$New_method/
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE+PCA/$New_method/
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE/$New_method/
mkdir MRNA_pipe/Pipe/rnaSeq/output/cluster/$New_method/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE+PCA/$New_method/$DATA_NAME
mkdir MRNA_pipe/Pipe/rnaSeq/output/visualize/tSNE/$New_method/$DATA_NAME