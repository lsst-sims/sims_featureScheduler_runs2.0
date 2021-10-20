#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 2
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -J presto_half
#SBATCH --ntasks-per-node=6

srun --output presto_half.out --ntasks=12 --multi-prog presto_cmds.conf
