#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 11
#SBATCH -t 24:00:00
#SBATCH -J presto_half

srun --output presto_half.out --ntasks=11 --multi-prog presto_cmds.conf
