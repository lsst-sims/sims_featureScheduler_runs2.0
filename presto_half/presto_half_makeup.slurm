#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH -J presto_half_makeup
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=20G   # memory per cpu-core

#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=yoachim@uw.edu

srun --output presto_half_makeup.out --ntasks=2 --multi-prog presto_makeup_cmds.conf
