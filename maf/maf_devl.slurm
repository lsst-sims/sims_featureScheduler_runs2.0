#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -J maf_all
#SBATCH --ntasks-per-node=20

cat maf_cmds.sh | parallel -j 20
