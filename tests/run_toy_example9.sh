#!/bin/bash

python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_ip_investigation --test_name toy_example9_100kworms_10 --n_worms 1E5 --cpus 10 --init_pop_strains 1E5 --mutation_rate 1E-4 --infection_sym_count 200 --growth_factor 1.5 --steady_state_runtime 1E5 --pop_size_thr 1E10 --escape_rate 1E-2
