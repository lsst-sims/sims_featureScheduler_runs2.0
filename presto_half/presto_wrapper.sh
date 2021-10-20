#!/bin/bash

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin

root_dir="/scratch/yoachim/git_repos/sims_featureScheduler_runs2.0/presto_half/"

python ${root_dir}presto_half.py --outDir ${root_dir} $1 
