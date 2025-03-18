#### Run the Wright-Fisher model at individual or population level ####

import numpy as np
from project_package.simplify import remove_empty_leaves_and_rescale_edges,remove_empty_leaves, subsample_pop, merge_graphs
from project_package.update_pop import update_pop,update_pop_with_tracers,update_pop2_with_tracers,update_pop2,update_pop3
from project_package.generate_pop import SymPop

import time
import multiprocessing

def run_generation_of_host_pop(freelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                     mutation_rate, steady_state_runtime,
                     max_runtime, growth_factor,pop_size_thr, 
                     stop_when_fixed=True, simplify=1,
                     verbose=0, t=0,sampling_rate=np.inf,nthreads=1):
    
    Graph_list_hostassociated=[]
    Graph_list_escapees=[]

    ### parallel phase ###
    
    args=[{'freelivingG':freelivingG,'host_pop_gen':host_pop_gen, 'host_id':host_id,
          'infection_sym_count':infection_sym_count, 'steady_state_runtime':steady_state_runtime,
          'mutation_rate':mutation_rate,'max_runtime':max_runtime, 
          'growth_factor':growth_factor,'pop_size_thr':pop_size_thr,
          'escape_rate':escape_rate,'sampling_rate':1,'verbose':verbose} for host_id in range(n_worms)]

    with multiprocessing.Pool(processes=nthreads) as pool:
        Graph_lists=pool.map(grow_and_steady_from_freeliving,args)
    
    Graph_list_hostassociated,Graph_list_escapees=map(list, zip(*Graph_lists)) 
    
    # get the host-associated symbiont meta-pop

    merged_Graph_hostassociated=merge_graphs(Graph_list_hostassociated)
    merged_Graph_hostassociated=remove_empty_leaves_and_rescale_edges(merged_Graph_hostassociated)

    # get the new free-living pop (after symbionts escape worm)
    merged_Graph_freeliving=merge_graphs([freelivingG]+Graph_list_escapees)
    merged_Graph_freeliving=remove_empty_leaves_and_rescale_edges(merged_Graph_freeliving)


    return(merged_Graph_hostassociated, merged_Graph_freeliving)


def grow_and_steady_from_freeliving(args):

    '''This function updates one symbiont population within one host with a growing phase (rate defined bygrowth_factor) until it reaches a population size threshold (pop_size_thr). Then, it updates the population under a steady-state model where the population size is stable for another n generations (steady_state_runtime). 
    The modifiable arguments are:
    host_id (int) : the id of the individual host,
    freelivingG (networkx/SymPop object) : the initial free-living population,
    infection_sym_count (int) : the number of bacterial cells that infect the host,
    host_pop_gen (int) : the generation of host population,
    escape_rate (float) : the proportion of bacterial cells that can escape the host,
    mutation_rate (float) : the mutaiton rate per bacterial cell per bacterial "generation", 
    steady_state_runtime (int) : the number of bacterial generations after the intra-host bacterial population reaches its maximum size,
    max_runtime (int) : the very maximum number of generations to run; this is mostly to avoid getting stuck, 
    growth_factor (float) : the growth factor for the symbiont population when it is growing; 1.05 meaning 5% growth at eaxh generation,
    pop_size_thr (int) : maximum symbiont population size in the host,
    verbose (int) : frequency at which population information are printed,
    sampling_rate (int, default=1) : rate at which intra-host bacterial populations are sampled
    
    The fixed arguments are:
    stop_when_fixed=False : continue the population updating process even when the population is fixed,
    t=0 : index of the initial bacterial generation within the host
    simplify=1 : the level of simplification applied to the population graph after earch symbiont generation; 1 means the empty leaves are removed'''

    # Load arguments
    np.random.seed()
    
    host_id=args['host_id']
    freelivingG=args['freelivingG']
    infection_sym_count=args['infection_sym_count']
    host_pop_gen=args['host_pop_gen']
    escape_rate=args['escape_rate']
    mutation_rate=args['mutation_rate']
    steady_state_runtime=args['steady_state_runtime']
    max_runtime=args['max_runtime']
    growth_factor=args['growth_factor']
    pop_size_thr=args['pop_size_thr']
    verbose=args['verbose']
    sampling_rate=args['sampling_rate']
    t=0
    simplify=1

    # Infection: set new host and subsample bacteria from free-living population
    subsampleG=subsample_pop(freelivingG,infection_sym_count)
    new_avail_id=str(host_pop_gen)+'.'+str(host_id)+'.0'

    # Grow host
    results=run_until_fixation3(subsampleG, mutation_rate, max_runtime, new_avail_id, growth_factor=growth_factor, stop_when_fixed=False, pop_size_thr=pop_size_thr, simplify=1, verbose=verbose, t=0,sampling_rate=sampling_rate)
    
    # Keep steady population in host for another n (steady_state_runtime) generations
    t=list(results.keys())[-1]
    
    results.update(run_until_fixation3(results[list(results.keys())[-1]], mutation_rate, t+steady_state_runtime, new_avail_id, growth_factor=1, stop_when_fixed=False, pop_size_thr=np.inf,simplify=1, verbose=verbose, t=t, sampling_rate=sampling_rate))  

    # get summary results
    max_gen=len(results.keys())-1
    final_intrahost_pop=results[max_gen]
    final_escapees_pop=subsample_pop(results[max_gen],SymPop(results[max_gen]).pop_size * escape_rate)

    return(final_intrahost_pop,final_escapees_pop)


def grow_and_steady(subsampleG, mutation_rate, steady_state_runtime, max_runtime, new_avail_id, growth_factor, pop_size_thr, stop_when_fixed=True, simplify=1, verbose=0, t=0,sampling_rate=1):

    '''This function updates one symbiont population within one host with a growing phase (rate defined bygrowth_factor) until it reaches a population size threshold (pop_size_thr). Then, it updates the population under a steady-state model where the population size is stable for another n generations (steady_state_runtime)'''
    
    
    # Grow host
    results=run_until_fixation3(subsampleG, mutation_rate, max_runtime, new_avail_id, growth_factor=growth_factor, stop_when_fixed=True, pop_size_thr=pop_size_thr, simplify=simplify, verbose=verbose, t=t,sampling_rate=sampling_rate)
    
    # Keep steady population in host for another n (steady_state_runtime) generations
    t=list(results.keys())[-1]
    
    results.update(run_until_fixation3(results[list(results.keys())[-1]], mutation_rate, t+steady_state_runtime, new_avail_id, growth_factor=1, stop_when_fixed=True, pop_size_thr=np.inf,simplify=simplify, verbose=verbose, t=t, sampling_rate=sampling_rate))  
    
    return(results)
    
def run_until_fixation(G, mutation_rate, runtime, growth_factor=1, stop_when_fixed=True, verbose=0,myupdate_pop_function=update_pop):
    
    
    results={0:G}
    t=0
    while True:
        t+=1
        G_plus1=myupdate_pop_function(G,mutation_rate,growth_factor=growth_factor)
        G_plus1=remove_empty_leaves_and_rescale_edges(G_plus1)

        G=G_plus1.copy()
        results.update({t:G})
        
        if verbose>0:
            if t%verbose==0:
                print(t)
        if stop_when_fixed==True:
            print('Population has fixed at t=',t,' generations.')
            if len((results[t]).nodes())==1:
                return(results)
                break
        
        if t>runtime:
            print('max runtime reached: ', runtime)
            return(results)
            break    
            
def run_until_fixation3(G, mutation_rate, runtime, new_avail_id, growth_factor=1, stop_when_fixed=False, pop_size_thr=np.inf, simplify=0, verbose=0, t=0, sampling_rate=50):

    ''' Run until fixation for update_pop3 code'''
    
    
    results={t:G}
    
    while True:
        t+=1
        (G_plus1,new_avail_id)=update_pop3(G,mutation_rate, new_avail_id, growth_factor=growth_factor)

        if simplify==1:
            G_plus1=remove_empty_leaves(G_plus1)

        if simplify==2:
            G_plus1=remove_empty_leaves_and_rescale_edges(G_plus1)

        G=G_plus1.copy()

        if t%sampling_rate==0:
            results.update({t:G})
        
        if verbose>0:
            if t%verbose==0:
                print(t)
        if stop_when_fixed==True:
            if len(G.nodes())==1:
                # print('Population has fixed at t=',t,' generations.')
                results.update({t:G})
                return(results)
                break

        if sum([attr['abundance'] for node,attr in G.nodes(data=True)]) > pop_size_thr:
            # print( 'max population size (',str(pop_size_thr),') reached at ',str(t),' generations')
            results.update({t:G})
            return(results)
            break
        
        if t>runtime:
            # print( 'max runtime reached at ', str(runtime), ' generations')
            results.update({t:G})
            return(results)
            break 

def run_until_fixation3_partitionned_trophosome(G, mutation_rate, runtime, new_avail_id, growth_factor=1, stop_when_fixed=True, verbose=0):
    ''' Run until fixation for update_pop3 code. Only removes leaf nodes that have abundances=0 at each generation. This is so that graphs from parallel runs are easier to merge'''
    
    results={0:G}
    t=0
    while True:
        t+=1
        (G_plus1,new_avail_id)=update_pop3(G,mutation_rate, new_avail_id, growth_factor=growth_factor)
        G_plus1=remove_empty_leaves_and_rescale_edges(G_plus1)

        G=G_plus1.copy()
        results.update({t:G})
        
        if verbose>0:
            if t%verbose==0:
                print(t)
        if stop_when_fixed==True:
            if len((results[t]).nodes())==1:
                return(results)
                break
        
        if t>runtime:
            return(results)
            break    