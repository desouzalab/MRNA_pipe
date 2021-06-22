#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=200000M

module load r/4.0.0

Rscript ../../rnaSeq/code/cluster/${CLUSTERMETHOD}/${CLUSTERMETHOD}.R --data_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --console_output_directory ../../rnaSeq/output/cluster/${CLUSTERMETHOD}/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --name_dataset ${DATASET}

