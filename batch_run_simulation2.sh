#!/bin/bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH -p public
#SBATCH -q public
#SBATCH -e slurm/slurm.%j.err
#SBATCH -o slurm/slurm.%j.out
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load mamba/latest
source activate python3_10

python ~/PyramidalCellSimulations/batch_run_simulation2.py
