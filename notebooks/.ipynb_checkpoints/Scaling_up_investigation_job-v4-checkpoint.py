#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import pickle
import pandas as pd
import sys
import os
# sys.path.append("/Users/maeva/Desktop/Modeling_trophosome/src/")
sys.path.append("/home/qiulab/data/CRF_project/work/Modeling_trophosome/src/")

from project_package.generate_pop import generate_initial_pop_unlinked, generate_random_fisherlog_pop_unlinked,generate_random_fisherlog_pop_binomial_tree, SymPop
from project_package.update_pop import update_pop3
from project_package.run_model import run_generation_of_host_pop
from project_package.plot import visualize_pop
from project_package.simplify import merge_graphs

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from alive_progress import alive_bar
import datetime
import time
import warnings
####################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Scaling_up_investigation"
    )
    parser.add_argument("--test_name", required=True, type=str)
    parser.add_argument("--mutation_rate", required=False, default=1E-12, type=lambda x: float(x))
    parser.add_argument("--growth_factor", required=False, default=1.2, type=lambda x: float(x))
    parser.add_argument("--steady_state_runtime", required=False, default=50, type=lambda x: int(float(x)))
    parser.add_argument("--max_runtime", required=False, default=np.inf)
    parser.add_argument("--pop_size_thr", required=False, default=1E4, type=lambda x: int(float(x)))
    parser.add_argument("--simplify", required=False, default=1, type=int)
    parser.add_argument("--verbose", required=False, default=0, type=int)
    parser.add_argument("--sampling_rate", required=False, default=1, type=int)
    parser.add_argument("--n_worms", required=True, default=10, type=lambda x: int(float(x)))
    parser.add_argument("--infection_sym_count", required=False, default=10, type=lambda x: int(float(x)))
    parser.add_argument("--tot_host_pop_gen", required=False, default=20, type=lambda x: int(float(x)))
    parser.add_argument("--escape_rate", required=False, default=1E-2, type=lambda x: float(x))
    parser.add_argument("--n_trials", required=False, default=1, type=int)
    parser.add_argument("--cpus", required=True, default=1, type=int)
    parser.add_argument("--init_pop_strains", required=False, default=100, type=lambda x: int(float(x)))
    parser.add_argument("--is_init_pop_log_fisher", required=False, default=True, type=bool)
    parser.add_argument("--out_dir", required=True, type=str)
    parser.add_argument("--random_seed", required=False, default=False, type=str)



    args = parser.parse_args()
    
    test=args.test_name
    out_dir=args.out_dir
    random_seed=args.random_seed
    
    ## init results objects    
    time_series_freeliving={}
    time_series_hostassociated={}
    
    ## params grow_and_steady
    mutation_rate=args.mutation_rate
    growth_factor=args.growth_factor
    steady_state_runtime=args.steady_state_runtime
    max_runtime=args.max_runtime
    pop_size_thr=args.pop_size_thr
    simplify=args.simplify
    verbose=args.verbose
    sampling_rate=args.sampling_rate
    
    ## params run_host_pop_gen
    n_worms=args.n_worms
    infection_sym_count=args.infection_sym_count
    tot_host_pop_gen=args.tot_host_pop_gen
    escape_rate=args.escape_rate
    n_trials=args.n_trials
    cpus=args.cpus
    
    # init tree log-fisher
    init_pop_strains=args.init_pop_strains
    is_init_pop_log_fisher=args.is_init_pop_log_fisher

####################

### print params   
    
myparams=['## params grow_and_steady',
'mutation_rate='+str(mutation_rate),
'growth_factor='+str(growth_factor),
'steady_state_runtime='+str(steady_state_runtime),
'max_runtime='+str(max_runtime),
'pop_size_thr='+str(pop_size_thr),
'simplify='+str(simplify),
'verbose='+str(verbose),
'sampling_rate='+str(sampling_rate),
'',
'## params run_host_pop_gen',
'n_worms='+str(n_worms),
'infection_sym_count='+str(infection_sym_count),
'tot_host_pop_gen='+str(tot_host_pop_gen),
'escape_rate='+str(escape_rate),
'n_trials='+str(n_trials),
'cpus='+str(cpus),
'',
'## init tree',
'init_pop_strains='+str(init_pop_strains),
'is_init_pop_log_fisher='+str(is_init_pop_log_fisher)]
    
dump_dir=out_dir+'/'+datetime.datetime.now().strftime('%Y_%m_%d')
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)
if os.path.isdir(dump_dir)==False:
    os.mkdir(dump_dir)
if os.path.isfile(dump_dir+'/'+'times.txt')==False:
    with open(dump_dir+'/'+'times.txt','w') as fa:
        fa.write('\t'.join(['test','n_worms','cpus','write_time','run_time'])+'\n')
    
with open(dump_dir+'/'+test+'.input_params.txt', 'w') as fa:
    fa.write('\n'.join(myparams))

##################

if random_seed==False:
    np.random.seed(666)

def fxn():
    warnings.warn("deprecated", DeprecationWarning)
    
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

        
    if is_init_pop_log_fisher==True:  
        freelivingG=generate_random_fisherlog_pop_unlinked(i=init_pop_strains) #n=individuals, i=strains
        freelivingG=nx.relabel_nodes(freelivingG, dict([[i,'0.0.'+str(i)] for i in range(init_pop_strains)]))
        print("Initial symbiont population follows log-fisher distribution. Number of strains = ",init_pop_strains," Number of individuals = ",SymPop(freelivingG).pop_size)
        
    ## init tree
    
    # fitnesses=[0.95,0.9] # allele_id:fitness
    # abundances=[300000000,700000000] # allele_id:abundance # total 1G cells
    # freelivingG=generate_initial_pop_unlinked(fitnesses,abundances)
    # freelivingG=nx.relabel_nodes(freelivingG, {0:'0.0.0',1:'0.0.1'})
    # time_series_freeliving[0]=freelivingG
    
        
    ###############################
    ## init results objects    
    time_series_freeliving={}
    time_series_hostassociated={}
    print('test = ',test,'n_worms = ',n_worms,'cpus = ',cpus)
    start_time=time.time()
    
    for trial in range(n_trials):
        host_pop_gen=0
        time_series_freeliving[trial]=[freelivingG]
        time_series_hostassociated[trial]=[nx.Graph()]
        myfreelivingG=freelivingG.copy()
    
        for host_pop_gen in range(1,tot_host_pop_gen+1):
        
            merged_Graph_hostassociated, merged_Graph_freeliving,_=run_generation_of_host_pop(myfreelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                                 mutation_rate, steady_state_runtime,
                                 max_runtime, growth_factor=growth_factor,
                                 stop_when_fixed=True, pop_size_thr=pop_size_thr, simplify=simplify,
                                 verbose=verbose, t=0,sampling_rate=sampling_rate,nthreads=cpus)
            
            time_series_hostassociated[trial]+=[merged_Graph_hostassociated]
            time_series_freeliving[trial]+=[merged_Graph_freeliving]
            myfreelivingG=merged_Graph_freeliving
                
    time_run=time.time()-start_time
    # time_df.loc[time_df.test==test,'time_run']=time_run
    ###########
    
    
    with open(dump_dir+'/'+test+'.time_series_hostassociated.graphs.pickledump', 'wb') as handle:
      pickle.dump(time_series_hostassociated, handle)
    
    with open(dump_dir+'/'+test+'.time_series_freeliving.graphs.pickledump', 'wb') as handle:
      pickle.dump(time_series_freeliving, handle)

    time_write=time.time()

    appended_line='\t'.join([test,str(n_worms),str(cpus),str(time_write-time_run),str(time_run)])

    with open(dump_dir+'/'+'times.txt','a') as fa:
        fa.write(appended_line+'\n')



