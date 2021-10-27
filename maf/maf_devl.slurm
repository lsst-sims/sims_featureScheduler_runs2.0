#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -J maf_all
#SBATCH --ntasks-per-node=20

srun --output maf_all.out --ntasks=127 --multi-prog maf_cmds.conf
