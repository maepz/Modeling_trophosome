#### Run the Wright-Fisher model at individual or population level ####

import numpy as np
from project_package.simplify import remove_empty_leaves_and_rescale_edges,remove_empty_leaves, subsample_pop,merge_graphs
from project_package.update_pop import update_pop,update_pop_with_tracers,update_pop2_with_tracers,update_pop2,update_pop3
from project_package.tskit_functions import prune_unique_genotypes_add_metadata,merge_by_root_and_combine_ancestrals,NodeTable_to_pandasDF,cell_coalescent_tree_with_mutations,tstree_to_networkx_graph
from project_package.generate_pop import SymPop
import networkx as nx
import time
import multiprocessing

import matplotlib.pyplot as plt
from project_package.plot import visualize_pop

def run_generation_of_host_pop_v4(freelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                     mutation_rate, steady_state_runtime,
                     max_runtime, growth_factor,pop_size_thr, 
                     stop_when_fixed=True, simplify=1,
                     verbose=0, t=0,sampling_rate=np.inf,nthreads=1):

    ''' Same as run_generation_of_host_pop but intra-host selection has been turned off
    '''
    ### Get free-living population paramaters

    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in freelivingG.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    
    ### parallel phase ###
    
    args=[{'freelivingG':freelivingG,'alleles':alleles,'abundances':abundances,'fitnesses':fitnesses,
           'host_pop_gen':host_pop_gen, 'host_id':host_id,
          'infection_sym_count':infection_sym_count, 'steady_state_runtime':steady_state_runtime,
          'mutation_rate':mutation_rate,'max_runtime':max_runtime, 
          'growth_factor':growth_factor,'pop_size_thr':pop_size_thr,
          'escape_rate':escape_rate,'sampling_rate':1,'verbose':verbose} for host_id in range(n_worms)]

    with multiprocessing.Pool(processes=nthreads) as pool:
        # print(args)
        Graph_lists=pool.map(grow_and_steady_from_freeliving_v4,args)
    
    Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees,results=map(list, zip(*Graph_lists)) 

    # get the initial host-associated symbiont meta-pop

    merged_Graph_list_init_intrahost=merge_graphs(Graph_list_init_intrahost)
    
    # get the host-associated symbiont meta-pop

    merged_Graph_hostassociated=merge_graphs(Graph_list_hostassociated)
    merged_Graph_hostassociated=remove_empty_leaves_and_rescale_edges(merged_Graph_hostassociated)

    # get the new free-living pop (after symbionts escape worm)
    init_pop_attr=np.array([[node,attr['abundance']]for node,attr in merged_Graph_list_init_intrahost.nodes(data=True) if attr['abundance']>0])
    init_pop_alleles = init_pop_attr[:,0]
    init_pop_abundances = np.array(list(map(int,init_pop_attr[:,1])))
        
    adj_freelivingG=[[alleles[i],{'abundance':abundances[i]-init_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(init_pop_abundances)) if init_pop_abundances[i] > 0]

    freelivingG.update(nodes=adj_freelivingG)
    
    merged_Graph_freeliving=merge_graphs([freelivingG]+Graph_list_escapees)
    merged_Graph_freeliving=remove_empty_leaves_and_rescale_edges(merged_Graph_freeliving)

    return(merged_Graph_hostassociated, merged_Graph_freeliving,results)

def run_generation_of_host_pop_v3(freelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                     mutation_rate, pop_size_thr, 
                     t=0,nthreads=1):
    
    '''this is the same as run_generation_of_host_pop but uses tskit instead of the Wright-Fisher to get Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees'''

    
    ### Get free-living population paramaters

    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in freelivingG.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    
    ### parallel phase ###
    
    args=[{'freelivingG':freelivingG,'alleles':alleles,'abundances':abundances,'fitnesses':fitnesses,
           'host_pop_gen':host_pop_gen, 'host_id':host_id,
          'infection_sym_count':infection_sym_count, 
          'mutation_rate':mutation_rate,'pop_size_thr':pop_size_thr,
          'escape_rate':escape_rate} for host_id in range(n_worms)]


    with multiprocessing.Pool(processes=nthreads) as pool:
        Graph_lists=pool.map(grow_from_freeliving,args)
        
    Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees=map(list, zip(*Graph_lists)) 


    # #### looped ###
    # Graph_list_init_intrahost=[]
    # Graph_list_hostassociated=[]
    # Graph_list_escapees=[]
    # for arg in args:
    #     print('host_id',arg['host_id'])
    #     res1,res2,res3=grow_from_freeliving(arg)
    #     Graph_list_init_intrahost+=[res1]
    #     Graph_list_hostassociated+=[res2]
    #     Graph_list_escapees+=[res3]

    
    # get the initial host-associated symbiont meta-pop

    merged_Graph_list_init_intrahost=merge_graphs(Graph_list_init_intrahost)
    
    # get the host-associated symbiont meta-pop

    merged_Graph_hostassociated=merge_graphs(Graph_list_hostassociated)
    merged_Graph_hostassociated=remove_empty_leaves_and_rescale_edges(merged_Graph_hostassociated)

    # get the new free-living pop (after symbionts escape worm)
    init_pop_attr=np.array([[node,attr['abundance']]for node,attr in merged_Graph_list_init_intrahost.nodes(data=True) if attr['abundance']>0])
    init_pop_alleles = init_pop_attr[:,0]
    init_pop_abundances = np.array(list(map(int,init_pop_attr[:,1])))
        
    adj_freelivingG=[[alleles[i],{'abundance':abundances[i]-init_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(init_pop_abundances)) if init_pop_abundances[i] > 0]

    freelivingG.update(nodes=adj_freelivingG)
    
    merged_Graph_freeliving=merge_graphs([freelivingG]+Graph_list_escapees)
    merged_Graph_freeliving=remove_empty_leaves_and_rescale_edges(merged_Graph_freeliving)

    return(merged_Graph_hostassociated, merged_Graph_freeliving)
    
def run_generation_of_host_pop_v2(freelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                     mutation_rate, pop_size_thr, 
                     t=0,nthreads=1):
    
    '''This is the same as run_generation_of_host_pop but uses tskit instead of the Wright-Fisher to get Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees. Ploidy=1 and growth rate are set.'''

    
    ### Get free-living population paramaters

    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in freelivingG.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    
    ### parallel phase ###
    
    args=[{'freelivingG':freelivingG,'alleles':alleles,'abundances':abundances,'fitnesses':fitnesses,
           'host_pop_gen':host_pop_gen, 'host_id':host_id,
          'infection_sym_count':infection_sym_count, 
          'mutation_rate':mutation_rate,'pop_size_thr':pop_size_thr,
          'escape_rate':escape_rate} for host_id in range(n_worms)]


    # with multiprocessing.Pool(processes=nthreads) as pool:
    #     Graph_lists=pool.map(grow_from_freeliving,args)
        
    # Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees=map(list, zip(*Graph_lists)) 


    #### looped ###
    Graph_list_init_intrahost=[]
    Graph_list_hostassociated=[]
    Graph_list_escapees=[]
    for arg in args:
        print('host_id',arg['host_id'])
        res1,res2,res3=grow_from_freeliving(arg)
        Graph_list_init_intrahost+=[res1]
        Graph_list_hostassociated+=[res2]
        Graph_list_escapees+=[res3]

    
    # get the initial host-associated symbiont meta-pop

    merged_Graph_list_init_intrahost=merge_graphs(Graph_list_init_intrahost)
    
    # get the host-associated symbiont meta-pop

    merged_Graph_hostassociated=merge_graphs(Graph_list_hostassociated)
    merged_Graph_hostassociated=remove_empty_leaves_and_rescale_edges(merged_Graph_hostassociated)

    # get the new free-living pop (after symbionts escape worm)
    init_pop_attr=np.array([[node,attr['abundance']]for node,attr in merged_Graph_list_init_intrahost.nodes(data=True) if attr['abundance']>0])
    init_pop_alleles = init_pop_attr[:,0]
    init_pop_abundances = np.array(list(map(int,init_pop_attr[:,1])))
        
    adj_freelivingG=[[alleles[i],{'abundance':abundances[i]-init_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(init_pop_abundances)) if init_pop_abundances[i] > 0]

    freelivingG.update(nodes=adj_freelivingG)
    
    merged_Graph_freeliving=merge_graphs([freelivingG]+Graph_list_escapees)
    merged_Graph_freeliving=remove_empty_leaves_and_rescale_edges(merged_Graph_freeliving)

    return(merged_Graph_hostassociated, merged_Graph_freeliving)

def run_generation_of_host_pop(freelivingG, n_worms, infection_sym_count,host_pop_gen,escape_rate,
                     mutation_rate, steady_state_runtime,
                     max_runtime, growth_factor,pop_size_thr, 
                     stop_when_fixed=True, simplify=1,
                     verbose=0, t=0,sampling_rate=np.inf,nthreads=1):

    # Graph_list_init_intrahost=[]
    # Graph_list_hostassociated=[]
    # Graph_list_escapees=[]
    
    ### Get free-living population paramaters

    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in freelivingG.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    
    ### parallel phase ###
    
    args=[{'freelivingG':freelivingG,'alleles':alleles,'abundances':abundances,'fitnesses':fitnesses,
           'host_pop_gen':host_pop_gen, 'host_id':host_id,
          'infection_sym_count':infection_sym_count, 'steady_state_runtime':steady_state_runtime,
          'mutation_rate':mutation_rate,'max_runtime':max_runtime, 
          'growth_factor':growth_factor,'pop_size_thr':pop_size_thr,
          'escape_rate':escape_rate,'sampling_rate':1,'verbose':verbose} for host_id in range(n_worms)]

    with multiprocessing.Pool(processes=nthreads) as pool:
        # print(args)
        Graph_lists=pool.map(grow_and_steady_from_freeliving,args)
    
    Graph_list_init_intrahost,Graph_list_hostassociated,Graph_list_escapees,results=map(list, zip(*Graph_lists)) 
    # print('len grapgs',len(Graph_list_init_intrahost), len(Graph_list_hostassociated), len(Graph_list_escapees))

    # get the initial host-associated symbiont meta-pop

    merged_Graph_list_init_intrahost=merge_graphs(Graph_list_init_intrahost)
    
    # get the host-associated symbiont meta-pop

    merged_Graph_hostassociated=merge_graphs(Graph_list_hostassociated)
    merged_Graph_hostassociated=remove_empty_leaves_and_rescale_edges(merged_Graph_hostassociated)

    # get the new free-living pop (after symbionts escape worm)
    init_pop_attr=np.array([[node,attr['abundance']]for node,attr in merged_Graph_list_init_intrahost.nodes(data=True) if attr['abundance']>0])
    init_pop_alleles = init_pop_attr[:,0]
    init_pop_abundances = np.array(list(map(int,init_pop_attr[:,1])))
    
    # print('init_pop [allele, abundance]',list(zip(init_pop_alleles,init_pop_abundances)))
    
    adj_freelivingG=[[alleles[i],{'abundance':abundances[i]-init_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(init_pop_abundances)) if init_pop_abundances[i] > 0]

    # print('freelivingG',freelivingG.nodes(data=True))
    # print('adj_freelivingG',adj_freelivingG)
    freelivingG.update(nodes=adj_freelivingG)
    # print('updated_freelivingG',freelivingG.nodes(data=True))
    # print('Graph_list_escapees',Graph_list_escapees)
    # print('Graph_list_escapees',[Graph_list_escapees[i].nodes(data=True) for i in range(len(Graph_list_escapees))])
    
    merged_Graph_freeliving=merge_graphs([freelivingG]+Graph_list_escapees)
    merged_Graph_freeliving=remove_empty_leaves_and_rescale_edges(merged_Graph_freeliving)
    # print('gen1_freelivingG',merged_Graph_freeliving.nodes(data=True))


    return(merged_Graph_hostassociated, merged_Graph_freeliving,results)

def grow_from_freeliving(args):

    '''This function grows one symbiont population within one host until it reaches a population size threshold (pop_size_thr). The cell
    phylogeny is reconstructed backwards through coalescence using tskit tree structure. Mutations are added to the populaiton tree.
    The modifiable arguments are:
    host_id (int) : the id of the individual host,
    alleles (list of str): the alleles of the free-living population
    abundances (list of int): the abundances of the free-living population
    fitnesses (list of floats): the fitnesses of the free-living population
    freelivingG (networkx/SymPop object) : the initial free-living population,
    infection_sym_count (int) : the number of bacterial cells that infect the host,
    host_pop_gen (int) : the generation of host population,
    escape_rate (float) : the proportion of bacterial cells that can escape the host,
    mutation_rate (float) : the mutaiton rate per bacterial cell per bacterial "generation", 
    pop_size_thr (int) : maximum symbiont population size in the host,
    
    '''
    
    # Load arguments
    # np.random.seed()
    
    host_id=args['host_id']
    freelivingG=args['freelivingG']
    alleles=args['alleles']
    abundances=args['abundances']
    fitnesses=args['fitnesses']
    infection_sym_count=args['infection_sym_count']
    host_pop_gen=args['host_pop_gen']
    escape_rate=args['escape_rate']
    mutation_rate=args['mutation_rate']
    pop_size_thr=args['pop_size_thr']
    seed=None

    
    tree_seed = np.random.seed(seed)
    mutation_seed = np.random.seed(seed)
    fitness_seed = np.random.seed(seed)

    # Infection: set new host and subsample bacteria from free-living population
    
    initial_intrahost_pop=subsample_pop(freelivingG,infection_sym_count)
    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in initial_intrahost_pop.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    
    new_avail_id=str(host_pop_gen)+'.'+str(host_id)+'.0'

    # Grow host: V1: all cells that enter the host ends up in the adult in the same proportions. Will have to play with msprime.sim_ancestry for different demographic models. https://tskit.dev/msprime/docs/stable/api.html#msprime.sim_ancestry
    
    G_list=[]
    for i in range(len(alleles)):
        allele=alleles[i]
        abundance=abundances[i]
        fitness=fitnesses[i]
    
        #For each cell that infected the host, trace back tree and mutations, put the trees into an array
        counter=1
        mts_unique_list=[]
        while counter <= abundance: #coalesce:
            n=pop_size_thr/infection_sym_count   # here we set the size of the cell line population in the final adult    
            mts_unique=cell_coalescent_tree_with_mutations(n=n, mu=mutation_rate, fitness=fitness, tree_seed=tree_seed, mutation_seed=mutation_seed, fitness_seed=fitness_seed)
            mts_unique_list+=[mts_unique]
            if counter == abundance:
                strain_tree=merge_by_root_and_combine_ancestrals(mts_unique_list)
                G=tstree_to_networkx_graph(strain_tree) 
        
                rootnode=[i for i in range(len(strain_tree.tables.nodes)) if (strain_tree.node(i).metadata.get("abundance", 0) > 0) and (strain_tree.node(i).metadata.get("mutation_count", 0) == 0) ]
                new_avail_prefix='.'.join(new_avail_id.split('.')[0:-1])
                new_avail_idx=int(new_avail_id.split('.')[-1])
                new_alleles_ids=['.'.join([new_avail_prefix,str(i)]) for i in range(new_avail_idx,new_avail_idx+len(G.nodes))]
                new_avail_id='.'.join([new_avail_prefix,str(new_avail_idx+len(G.nodes)-1)])
                
                node_id_dic=dict(zip(G.nodes,new_alleles_ids))
                node_id_dic[rootnode[0]]=str(allele)
                
                nx.relabel_nodes(G, node_id_dic, copy=False)
                G_list+=[G]
                
            counter+=1
    
    
            
    final_intrahost_pop=nx.union_all(G_list)   
    final_escapees_pop=subsample_pop(final_intrahost_pop,SymPop(final_intrahost_pop).pop_size * escape_rate)
    
    return(initial_intrahost_pop,final_intrahost_pop,final_escapees_pop)

def grow_and_steady_from_freeliving_v4(args):

    '''This function updates one symbiont population within one host with a growing phase (rate defined bygrowth_factor) until it reaches a population size threshold (pop_size_thr). Then, it updates the population under a steady-state model where the population size is stable for another n generations (steady_state_runtime). As defined by run_until_fixation4 (running update_pop4); there is NO selection within the host.
    The modifiable arguments are:
    host_id (int) : the id of the individual host,
    alleles (list of str): the alleles of the free-living population
    abundances (list of int): the abundances of the free-living population
    fitnesses (list of floats): the fitnesses of the free-living population
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
    alleles=args['alleles']
    abundances=args['abundances']
    fitnesses=args['fitnesses']
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
    #subsampleG=subsample_pop(freelivingG,infection_sym_count)
    
    tot=sum(abundances)
    weights=np.multiply(abundances,1/tot)
    new_pop_abundances = np.random.multinomial(int(infection_sym_count), weights)
    # print(new_pop_abundances)
    new_pop_alleles = [alleles[i] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    
    subsampleG=nx.Graph(freelivingG.subgraph(new_pop_alleles))
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    subsampleG.update(nodes=adj)
    initial_intrahost_pop=subsampleG.copy()
    
    new_avail_id=str(host_pop_gen)+'.'+str(host_id)+'.0'

    # Grow host
    results=run_until_fixation4(subsampleG, mutation_rate, max_runtime, new_avail_id, growth_factor=growth_factor, stop_when_fixed=False, pop_size_thr=pop_size_thr, simplify=1, verbose=verbose, t=0,sampling_rate=sampling_rate)
    
    # Keep steady population in host for another n (steady_state_runtime) generations
    t=list(results.keys())[-1]
    # print('t',t)
    
    results.update(run_until_fixation4(results[list(results.keys())[-1]], mutation_rate, t+steady_state_runtime, new_avail_id, growth_factor=1, stop_when_fixed=False, pop_size_thr=np.inf,simplify=1, verbose=verbose, t=t, sampling_rate=sampling_rate))  

    # get summary results
    max_gen=len(results.keys())-1
    final_intrahost_pop=results[max_gen]
    # print('SymPop(results[max_gen]).pop_size',SymPop(results[max_gen]).pop_size)
    # print('escapee_count',SymPop(results[max_gen]).pop_size * escape_rate)
    final_escapees_pop=subsample_pop(final_intrahost_pop,SymPop(results[max_gen]).pop_size * escape_rate)
    return(initial_intrahost_pop,final_intrahost_pop,final_escapees_pop,results)
    
def grow_and_steady_from_freeliving(args):

    '''This function updates one symbiont population within one host with a growing phase (rate defined bygrowth_factor) until it reaches a population size threshold (pop_size_thr). Then, it updates the population under a steady-state model where the population size is stable for another n generations (steady_state_runtime). As defined by run_until_fixation3 (running update_pop3), the symbionts passing into the next generation are sampled using a multinomial which is weighed by fitness (i.e. there is selection within the host).
    The modifiable arguments are:
    host_id (int) : the id of the individual host,
    alleles (list of str): the alleles of the free-living population
    abundances (list of int): the abundances of the free-living population
    fitnesses (list of floats): the fitnesses of the free-living population
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
    alleles=args['alleles']
    abundances=args['abundances']
    fitnesses=args['fitnesses']
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
    #subsampleG=subsample_pop(freelivingG,infection_sym_count)
    
    tot=sum(abundances)
    weights=np.multiply(abundances,1/tot)
    new_pop_abundances = np.random.multinomial(int(infection_sym_count), weights)
    # print(new_pop_abundances)
    new_pop_alleles = [alleles[i] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    
    subsampleG=nx.Graph(freelivingG.subgraph(new_pop_alleles))
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    subsampleG.update(nodes=adj)
    initial_intrahost_pop=subsampleG.copy()
    
    new_avail_id=str(host_pop_gen)+'.'+str(host_id)+'.0'

    # Grow host
    results=run_until_fixation3(subsampleG, mutation_rate, max_runtime, new_avail_id, growth_factor=growth_factor, stop_when_fixed=False, pop_size_thr=pop_size_thr, simplify=1, verbose=verbose, t=0,sampling_rate=sampling_rate)
    
    # Keep steady population in host for another n (steady_state_runtime) generations
    t=list(results.keys())[-1]
    # print('t',t)
    
    results.update(run_until_fixation3(results[list(results.keys())[-1]], mutation_rate, t+steady_state_runtime, new_avail_id, growth_factor=1, stop_when_fixed=False, pop_size_thr=np.inf,simplify=1, verbose=verbose, t=t, sampling_rate=sampling_rate))  

    # get summary results
    max_gen=len(results.keys())-1
    final_intrahost_pop=results[max_gen]
    # print('SymPop(results[max_gen]).pop_size',SymPop(results[max_gen]).pop_size)
    # print('escapee_count',SymPop(results[max_gen]).pop_size * escape_rate)
    final_escapees_pop=subsample_pop(final_intrahost_pop,SymPop(results[max_gen]).pop_size * escape_rate)
    return(initial_intrahost_pop,final_intrahost_pop,final_escapees_pop,results)


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

def run_until_fixation4(G, mutation_rate, runtime, new_avail_id, growth_factor=1, stop_when_fixed=False, pop_size_thr=np.inf, simplify=0, verbose=0, t=0, sampling_rate=50):

    ''' Run until fixation for update_pop4 code'''
    
    
    results={t:G}
    
    while True:
        t+=1
        (G_plus1,new_avail_id)=update_pop4(G,mutation_rate, new_avail_id, growth_factor=growth_factor)

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