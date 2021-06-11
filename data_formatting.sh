#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=200000M

conda activate env_full

Python data_fromatting.py