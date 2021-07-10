

#!/bin/bash
DATASET=GSE74672
sbatch --job-name=${DATASET} --output=/home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/cluster/backspin/${DATASET}/CALL_SPIN.out --export=DATASET=${DATASET} /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/cluster/backspin/call_spin.sh
