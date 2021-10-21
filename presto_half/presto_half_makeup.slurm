#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 12:00:00
#SBATCH -J presto_half_makeup
#SBATCH --ntasks-per-node=2

srun --output presto_half_makeup.out --ntasks=2 --multi-prog presto_makeup_mds.conf
