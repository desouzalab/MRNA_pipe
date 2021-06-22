#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8000M

module load r/4.0.0

Rscript ../../rnaSeq/code/visualize/${VISUALIZEMETHOD}/${VISUALIZEMETHOD}.R --output_directory ../../rnaSeq/output/visualize/${VISUALIZEMETHOD}/${CLUSTERMETHOD}/${DATASET}/ --preprocessed_input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}/ --trueCluster_input_directory ~/projects/def-cdesouza/Lab/data/true/${DATASET}/ --cluster_input_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET}/ --name_dataset ${DATASET}
