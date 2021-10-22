#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -J presto_half_maf
#SBATCH --ntasks-per-node=20

srun --output maf.out --ntasks=12 --multi-prog maf_cmds.conf
