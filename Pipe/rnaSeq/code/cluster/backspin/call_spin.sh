#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=00:60:00
#SBATCH --mem-per-cpu=200000M

backspin -i /home/emiliano/projects/def-cdesouza/Lab/data/preprocessed/${DATASET}/ceftype_1_${DATASET}_inputforBackSpin.cef -o /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/cluster/backspin/${DATASET}/backspin.cef -f 500 -v -d 6