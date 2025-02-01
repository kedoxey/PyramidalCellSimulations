#!/bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
#SBATCH -p general
#SBATCH -q public
#SBATCH -e logs/slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load mamba/latest
source activate python3_10

python ~/PyramidalCellSimulations/batch_run_simulation.py
