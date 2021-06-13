#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=200000M
#sbatch --parsable --job-name=pre --output=/home/emiliano/projects/def-cdesouza/Lab/data/pythonLog.out /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/data_formatting.sh
conda activate env_full

python data_fromatting.py