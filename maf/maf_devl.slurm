#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 127
#SBATCH -t 48:00:00
#SBATCH -J maf_all
#SBATCH --ntasks-per-node=18

srun --output presto_half.out --ntasks=127 --multi-prog maf_cmds.conf
