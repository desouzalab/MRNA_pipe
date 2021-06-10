#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=4000M

module load r/4.0.0

Rscript ../../rnaSeq/code/preprocess/preprocess.R --data_output_directory ~/projects/def-cdesouza/Lab/data/preprocessed/${DATASET} --console_output_directory ../../rnaSeq/output/preprocess/${DATASET} --input_directory ~/projects/def-cdesouza/Lab/data/raw/${DATASET}/ --name_dataset ${DATASET}

