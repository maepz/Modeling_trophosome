import numpy as np
import time

def update_pop3(G,mutation_rate, new_avail_id, growth_factor=1):
    ''' This function takes a population under the form of a networkx object and updates it to a new population following a Fisher-Wright process where the whole population is replaced by descendents and new alleles can arrise according to a specified mutation rate. This function is based un update_pop2 but add the seed/trial/host number to the strain name so two independent trials can be merged.
    OPTIONS:
    growth_factor=1 Growth factor for the population. Varies from 0 to Inf. [0-1] population reduces, [1-Inf] population expands.
    '''
    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in G.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))
    n=sum(abundances)

    abs_weights=np.multiply(abundances,fitnesses)
    sum_abs_wt=sum(abs_weights)
    rel_weights=np.multiply(abs_weights,1/sum_abs_wt)
    weights=rel_weights

    ########
    
    new_pop_abundances = np.random.multinomial(int(n * growth_factor), weights)
    tot_new_pop=sum(new_pop_abundances)
    mutated_ind_count=np.random.binomial(int(n), mutation_rate)

    if mutated_ind_count>0:
        mutated_ind = np.random.multinomial(mutated_ind_count,new_pop_abundances/tot_new_pop)
        new_pop_abundances=new_pop_abundances-mutated_ind
    
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(alleles))]
    G.update(nodes=adj)
            
    if mutated_ind_count>0: # if new alleles have arisen
            
        # Update pop with new alleles
        mutated_alleles=np.nonzero(mutated_ind)[0]
        mutated_alleles_count =len(mutated_alleles)

        new_avail_prefix='.'.join(new_avail_id.split('.')[0:-1])
        new_avail_idx=int(new_avail_id.split('.')[-1])
        new_alleles_ids=['.'.join([new_avail_prefix,str(i)]) for i in range(new_avail_idx,new_avail_idx+mutated_ind_count)]
        new_avail_id='.'.join([new_avail_prefix,str(new_avail_idx+mutated_ind_count)])
        
        new_alleles_parent_indices=np.repeat(mutated_alleles, mutated_ind[mutated_alleles]).astype(int)
        mutated_alleles_parent_fitness=fitnesses[new_alleles_parent_indices]
        new_alleles_fitnesses=np.sum([mutated_alleles_parent_fitness,np.random.normal(-0.01, 0.01,size=mutated_ind_count)], axis=0) # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01
        
        new_alleles_nodes=[[new_alleles_ids[i],{'abundance':1,'fitness':new_alleles_fitnesses[i]}] for i in range(len(new_alleles_fitnesses))]
        new_alleles_edges=np.array([alleles[new_alleles_parent_indices],new_alleles_ids,[{"distance":1}] * len(new_alleles_ids)]).transpose()

        G.update(edges=new_alleles_edges, nodes=new_alleles_nodes)
    
    return(G, new_avail_id)
    
def update_pop2(G,mutation_rate,growth_factor=1):
    ''' This function takes a population under the form of a networkx object and updates it to a new population following a Fisher-Wright process where the whole population is replaced by descendents and new alleles can arrise according to a specified mutation rate. This function is more memory efficient than update_pop
    OPTIONS:
    growth_factor=1 Growth factor for the population. Varies from 0 to Inf. [0-1] population reduces, [1-Inf] population expands.
    '''
    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in G.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = pop_attr[:,1]
    fitnesses = pop_attr[:,2]
    n=sum(abundances)

    abs_weights=np.multiply(abundances,fitnesses)
    sum_abs_wt=sum(abs_weights)
    rel_weights=np.multiply(abs_weights,1/sum_abs_wt)
    weights=rel_weights

    ########
    
    new_pop_abundances = np.random.multinomial(n*growth_factor, weights)
    tot_new_pop=sum(new_pop_abundances)
    mutated_ind_count=np.random.binomial(n, mutation_rate)

    if mutated_ind_count>0:
        mutated_ind = np.random.multinomial(mutated_ind_count,new_pop_abundances/tot_new_pop)
        new_pop_abundances=new_pop_abundances-mutated_ind
    
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(alleles))]
    G.update(nodes=adj)
            
    if mutated_ind_count>0: # if new alleles have arisen
            
        # Update pop with new alleles
            
        mutated_alleles=np.nonzero(mutated_ind)[0]
        mutated_alleles_count = len(mutated_alleles)

        lastid=int(max(alleles))
        new_alleles_ids=range(lastid+1,lastid+1+mutated_ind_count)
        new_alleles_parent_indices=np.repeat(mutated_alleles, mutated_ind[mutated_alleles]).astype(int)
        mutated_alleles_parent_fitness=fitnesses[new_alleles_parent_indices]
        new_alleles_fitnesses=np.sum([mutated_alleles_parent_fitness,np.random.normal(-0.01, 0.01,size=mutated_ind_count)], axis=0) # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01
        
        new_alleles_nodes=[[new_alleles_ids[i],{'abundance':1,'fitness':new_alleles_fitnesses[i]}] for i in range(len(new_alleles_fitnesses))]
        new_alleles_edges=np.array([alleles[new_alleles_parent_indices].astype(int),new_alleles_ids,[{"distance":1}] * len(new_alleles_ids)]).transpose()

        G.update(edges=new_alleles_edges, nodes=new_alleles_nodes)
    
    return(G)
    
def update_pop2_with_tracers(G,mutation_rate,growth_factor=1):
    ''' This function takes a population under the form of a networkx object and updates it to a new population following a Fisher-Wright process where the whole population is replaced by descendents and new alleles can arrise according to a specified mutation rate. This function includes time and memory tracers and is more computationally efficient than update_pop.
    OPTIONS:
    growth_factor=1 Growth factor for the population. Varies from 0 to Inf. [0-1] population reduces, [1-Inf] population expands.
    '''

    start_time=time.time() ######

    # newG=G.copy()

    time_newG_copy=time.time()-start_time #####
    
    ############### v2 ###############
    # population={node: attr['abundance'] for node,attr in G.nodes(data=True) if attr['abundance']>0}
    # fitness_dic = {node: attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0}

    # alleles = np.array(list(population.keys()))
    # abundances=np.array(list(population.values()))
    # fitnesses=np.array(list(fitness_dic.values()))
    # n=sum(population.values())

    # del population
    # del fitness_dic
    # alleles=np.array([node for node,attr in G.nodes(data=True) if attr['abundance']>0])
    # fitnesses=np.array([attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0]) 

    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in G.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = pop_attr[:,1]
    fitnesses = pop_attr[:,2]
    n=sum(abundances)

    time_pop_dic=time.time()-start_time-time_newG_copy

    abs_weights=np.multiply(abundances,fitnesses)
    sum_abs_wt=sum(abs_weights)
    rel_weights=np.multiply(abs_weights,1/sum_abs_wt)
    weights=rel_weights

    time_pop_weights=time.time()-start_time-time_newG_copy-time_pop_dic

    time_pop_attr=time.time()-start_time-time_newG_copy
    ########
    
    new_pop_abundances = np.random.multinomial(n*growth_factor, weights)
    tot_new_pop=sum(new_pop_abundances)
    mutated_ind_count=np.random.binomial(n, mutation_rate)

    time_new_pop_init=time.time()-start_time-time_newG_copy-time_pop_attr ####

    if mutated_ind_count>0:
        mutated_ind = np.random.multinomial(mutated_ind_count,new_pop_abundances/tot_new_pop)
        new_pop_abundances=new_pop_abundances-mutated_ind
    
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(alleles))]
    # newG.update(nodes=adj)
    G.update(nodes=adj)
    
    time_new_pop_adj=time.time()-start_time-time_newG_copy-time_pop_attr-time_new_pop_init #####
        
    if mutated_ind_count>0: # if new alleles have arisen
            
        # Update pop with new alleles
            
        mutated_alleles=np.nonzero(mutated_ind)[0]
        mutated_alleles_count = len(mutated_alleles)

        lastid=int(max(alleles))
        new_alleles_ids=range(lastid+1,lastid+1+mutated_ind_count)
        new_alleles_parent_indices=np.repeat(mutated_alleles, mutated_ind[mutated_alleles]).astype(int)
        # try:
        mutated_alleles_parent_fitness=fitnesses[new_alleles_parent_indices]
        # except TypeError:
        #     print(new_alleles_parent_indices)
        new_alleles_fitnesses=np.sum([mutated_alleles_parent_fitness,np.random.normal(-0.01, 0.01,size=mutated_ind_count)], axis=0) # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01
        
        new_alleles_nodes=[[new_alleles_ids[i],{'abundance':1,'fitness':new_alleles_fitnesses[i]}] for i in range(len(new_alleles_fitnesses))]
        # try:
        new_alleles_edges=np.array([alleles[new_alleles_parent_indices].astype(int),new_alleles_ids,[{"distance":1}] * len(new_alleles_ids)]).transpose()
        # except TypeError:
        #     print('new_alleles_parent_indices', new_alleles_parent_indices)
        #     print('alleles[new_alleles_parent_indices]', alleles[new_alleles_parent_indices])
            
        # newG.update(edges=new_alleles_edges, nodes=new_alleles_nodes)
        G.update(edges=new_alleles_edges, nodes=new_alleles_nodes)

    time_new_pop_mut=time.time()-start_time-time_newG_copy-time_pop_attr-time_new_pop_init-time_new_pop_adj #####
    
    # return(newG,time_newG_copy,time_pop_dic,time_pop_weights, time_pop_attr, time_new_pop_init, time_new_pop_adj, time_new_pop_mut)
    return(G,time_newG_copy,time_pop_dic,time_pop_weights, time_pop_attr, time_new_pop_init, time_new_pop_adj, time_new_pop_mut)


def update_pop_with_tracers(G,mutation_rate,growth_factor=1):
    ''' This function takes a population under the form of a networkx object and updates it to a new population following a Fisher-Wright process where the whole population is replaced by descendents and new alleles can arrise according to a specified mutation rate
    OPTIONS:
    growth_factor=1 Growth factor for the population. Varies from 0 to Inf. [0-1] population reduces, [1-Inf] population expands.
    '''

    start_time=time.time() ######

    newG=G.copy()

    time_newG_copy=time.time()-start_time #####

    ############### v1 #################
    
    fitness_dic={node: attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0}
    population={node: attr['abundance'] for node,attr in G.nodes(data=True) if attr['abundance']>0}
    alleles = list(population.keys())     
    n=sum(list(population.values()))

    time_pop_dic=time.time()-start_time-time_newG_copy
    
    abs_weights=[population[x]*fitness_dic[x] for x in alleles] # absolute fitness
    rel_weights=[wt/sum(abs_weights) for wt in abs_weights] # relative fitness
    weights=rel_weights
    
    time_pop_weights=time.time()-start_time-time_newG_copy-time_pop_dic

    time_pop_attr=time.time()-start_time-time_newG_copy
    ######### 
    
    new_pop = dict(zip(alleles,np.random.multinomial(n*growth_factor, weights)))
    tot_new_pop=sum(new_pop.values())
    mutated_ind_count=np.random.binomial(n, mutation_rate)

    time_new_pop_init=time.time()-start_time-time_newG_copy-time_pop_attr ####
    
    
    if mutated_ind_count>0:
        mutated_ind = dict(zip(new_pop.keys(),np.random.multinomial(mutated_ind_count, 
                                                                [v/tot_new_pop for v in new_pop.values()])))
        new_pop={k:new_pop[k]-mutated_ind[k] for k in new_pop.keys()}
        mutated_ind={k: mutated_ind[k] for k in mutated_ind.keys()}
        
    adj=[[node_id,{'abundance':val,'fitness':fitness_dic[node_id]}] for node_id,val in new_pop.items()]
    newG.update(nodes=adj)

    time_new_pop_adj=time.time()-start_time-time_newG_copy-time_pop_attr-time_new_pop_init #####
    
    if mutated_ind_count>0: # if new alleles have arisen
        
        # Update pop with new alleles
        new_alleles_ids=range(max(new_pop.keys())+1,max(new_pop.keys())+1+sum(mutated_ind.values()))
        new_alleles_fitnesses=[sum(x) for x in zip([fitness_dic[k] for k, v in mutated_ind.items() for _ in range(v)],
                                                   np.random.normal(-0.01, 0.01,size=sum(mutated_ind.values())))]
 # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01

        
        new_alleles_nodes=[[new_alleles_ids[i],{'abundance':1,'fitness':new_alleles_fitnesses[i]}] for i in range(len(new_alleles_fitnesses))]
        
        new_alleles_edges=list(zip([k for k, v in mutated_ind.items() for _ in range(v)],
                             new_alleles_ids,
                             [{"distance":1}] * len(new_alleles_ids)))

        newG.update(edges=new_alleles_edges, nodes=new_alleles_nodes)

    time_new_pop_mut=time.time()-start_time-time_newG_copy-time_pop_attr-time_new_pop_init-time_new_pop_adj #####
    
    return(newG,time_newG_copy,time_pop_dic,time_pop_weights, time_pop_attr, time_new_pop_init, time_new_pop_adj, time_new_pop_mut)
    
def update_pop(G,mutation_rate,growth_factor=1):
    ''' This function takes a population under the form of a networkx object and updates it to a new population following a Fisher-Wright process where the whole population is replaced by descendents and new alleles can arrise according to a specified mutation rate
    OPTIONS:
    growth_factor=1 Growth factor for the population. Varies from 0 to Inf. [0-1] population reduces, [1-Inf] population expands.
    '''
        
    newG=G.copy()
    
    fitness_dic={node: attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0}
    population={node: attr['abundance'] for node,attr in G.nodes(data=True) if attr['abundance']>0}
    
    n=sum(list(population.values()))

    alleles = list(population.keys())
    abs_weights=[population[x]*fitness_dic[x] for x in alleles] # absolute fitness
    rel_weights=[wt/sum(abs_weights) for wt in abs_weights] # relative fitness

    weights=rel_weights
    
    ######### 
    new_pop = dict(zip(alleles,np.random.multinomial(n*growth_factor, weights)))
    tot_new_pop=sum(new_pop.values())
    mutated_ind_count=np.random.binomial(n, mutation_rate)
    
    if mutated_ind_count>0:
        mutated_ind = dict(zip(new_pop.keys(),np.random.multinomial(mutated_ind_count, 
                                                                [v/tot_new_pop for v in new_pop.values()])))
        new_pop={k:new_pop[k]-mutated_ind[k] for k in new_pop.keys()}
        mutated_ind={k: mutated_ind[k] for k in mutated_ind.keys()}
        
    adj=[[node_id,{'abundance':val,'fitness':fitness_dic[node_id]}] for node_id,val in new_pop.items()]
    newG.update(nodes=adj)

    if mutated_ind_count>0: # if new alleles have arisen
        
        # Update pop with new alleles
        new_alleles_ids=range(max(new_pop.keys())+1,max(new_pop.keys())+1+sum(mutated_ind.values()))
        new_alleles_fitnesses=[sum(x) for x in zip([fitness_dic[k] for k, v in mutated_ind.items() for _ in range(v)],
                                                   np.random.normal(-0.01, 0.01,size=sum(mutated_ind.values())))]
 # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01

        
        new_alleles_nodes=[[new_alleles_ids[i],{'abundance':1,'fitness':new_alleles_fitnesses[i]}] for i in range(len(new_alleles_fitnesses))]
        
        new_alleles_edges=list(zip([k for k, v in mutated_ind.items() for _ in range(v)],
                             new_alleles_ids,
                             [{"distance":1}] * len(new_alleles_ids)))

        newG.update(edges=new_alleles_edges, nodes=new_alleles_nodes)

    return(newG)