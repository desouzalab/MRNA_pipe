

#!/bin/bash


DATASET=GSE74672
COMPAREMETHOD=AIR+VM+Purity

sbatch --job-name="${DATASET} ${COMPAREMETHOD}" --output=/home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/output/compare/${COMPAREMETHOD}/${DATASET}/ch_in.out  --export=DATASET=${DATASET},COMPAREMETHOD=${COMPAREMETHOD} /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code/compare/compare_job.sh