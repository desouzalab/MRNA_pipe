#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=20000M

module load r/4.0.0

Rscript ../../rnaSeq/code/compare/${COMPAREMETHOD}/${COMPAREMETHOD}.R --output_directory ../../rnaSeq/output/compare/${COMPAREMETHOD}/${DATASET} --preprocessed_input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}/ --true_cluster_input_directory ~/projects/def-cdesouza/Lab/data/true/${DATASET}/ --sc3_cluster_input_directory ../../rnaSeq/output/cluster/sc3/${DATASET}/ --seurat_cluster_input_directory ../../rnaSeq/output/cluster/seurat/${DATASET}/ --name_dataset ${DATASET}
