#!/bin/bash

source ~/.bashrc
mamba activate maeva_env

N=$1


### Toy_example9



#### toy_example1 to 8, 18: variation of mutation rate ####
#### {'1':1e-12,'2':1e-11,'3':1e-10,'4':1e-9,'5':1e-6,'6':1e-3,'7':1e-4,'8':1e-5, '18':1e-2} ####

#mu=$2


#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_ip_investigation --test_name V1_toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --mutation_rate ${mu}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_msprime.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --mutation_rate ${mu}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_msprime_v2_2.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --mutation_rate ${mu}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_v1_4.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation/V1_4 --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --mutation_rate ${mu} --steady_state_runtime 0 

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_v3_1.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation/V3_1 --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --mutation_rate ${mu} --steady_state_runtime 0 

####  toy_example10 to 17: variation of intra-host population ####
#### {'1':1e4,'10':1e5,'11':1e6,'12':1e7,'13':1e8,'14':1e9,'15':1e10,'16':1e11,'17':1e12} ####

  
intrahost_pop_size=$2

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name v1_toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --pop_size_thr ${intrahost_pop_size}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_msprime.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --pop_size_thr ${intrahost_pop_size}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_msprime_v2_2.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --pop_size_thr ${intrahost_pop_size}

#python /home/qiulab/data/CRF_project/work/Modeling_trophosome/notebooks/Scaling_up_investigation_job_v3_1.py --out_dir /home/qiulab/data/CRF_project/scratch/Modeling_trophosome/Scaling_up_investigation/V3_1 --test_name toy_example${N}_100worms_10cpus --n_worms 100 --cpus 10 --pop_size_thr ${intrahost_pop_size}