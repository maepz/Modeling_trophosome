#!/bin/bash

source ~/.bashrc
mamba activate maeva_env

N_WORMS=$1
CPUS=$2

python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example1_${N_WORMS}worms_${CPUS}cpus --n_worms $N_WORMS --cpus $CPUS
