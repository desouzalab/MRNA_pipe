#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=00:60:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=30G


module load r/4.0.0

Rscript ../../rnaSeq/code/compare/${COMPAREMETHOD}/${COMPAREMETHOD}.R --output_directory ../../rnaSeq/output/compare/${COMPAREMETHOD}/${DATASET} --preprocessed_input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}/ --true_cluster_input_directory ~/projects/def-cdesouza/Lab/data/true/${DATASET}/ --sc3_cluster_input_directory ../../rnaSeq/output/cluster/sc3/${DATASET}/ --seurat_cluster_input_directory ../../rnaSeq/output/cluster/seurat/${DATASET}/ --giniclust_cluster_input_directory ../../rnaSeq/output/cluster/giniclust/${DATASET}/ --backspin_cluster_input_directory /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/cluster/backspin/${DATASET} --name_dataset ${DATASET}
