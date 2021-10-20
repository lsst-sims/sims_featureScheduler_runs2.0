#!/bin/bash

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin

scimaf_dir --db $1
glance_dir --db $1
