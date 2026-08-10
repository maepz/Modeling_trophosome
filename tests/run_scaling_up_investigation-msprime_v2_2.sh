#!/bin/bash

source ~/.bashrc
mamba activate maeva_env

N_WORMS=$1
CPUS=$2

values=(10 10 10 4 10 1 10 2 100 40 100 10 100 10 100 20 100 2 100 4 1000 10 1000 40 1000 20 1000 10 1000 4 10000 20 10000 10)

for ((i=0; i<${#values[@]}; i+=2)); do
N_WORMS=${values[i]}
CPUS=${values[i+1]}
echo "NWORMS=$N_WORMS CPUS=$CPUS"

python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_msprime_v2_2.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example1_${N_WORMS}worms_${CPUS}cpus --n_worms $N_WORMS --cpus $CPUS

done