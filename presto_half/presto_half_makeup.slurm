#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J presto_half_makeup
#SBATCH --ntasks-per-node=1
#SBATCH --mem=35G   

#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=yoachim@uw.edu

srun --output presto_half_makeup.out --ntasks=1 --multi-prog presto_makeup_cmds.conf
