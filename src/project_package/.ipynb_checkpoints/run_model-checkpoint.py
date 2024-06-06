#### Run the Wright-Fisher model ####
import numpy as np
from project_package.simplify import remove_empty_leaves_and_rescale_edges


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

def split_pop_in_half(G):
    from collections import Counter
    import random
    
    fitness_dic={node: attr['fitness'] for node,attr in G.nodes(data=True)}
    parent_pop={node: attr['abundance'] for node,attr in G.nodes(data=True)}
    n=sum(parent_pop.values())

    rnd = Counter(random.sample(list(parent_pop.keys()), counts=map(int,list(parent_pop.values())), k=round(n/2)))
    new_pop1=Counter(parent_pop)
    new_pop1.subtract(rnd)
    
    new_pop2=Counter(parent_pop)
    new_pop2.subtract(new_pop1)

    
    newG1=G.copy()
    adj=[[node_id,{'abundance':val,'fitness':fitness_dic[node_id]}] for node_id,val in new_pop1.items()]
    newG1.update(nodes=adj)
    
    newG2=G.copy()
    adj=[[node_id,{'abundance':val,'fitness':fitness_dic[node_id]}] for node_id,val in new_pop2.items()]
    newG2.update(nodes=adj)

    return(newG1,newG2)


def run_until_fixation(G, mutation_rate, runtime, growth_factor=1, stop_when_fixed=True, verbose=0):
    
    
    results={0:G}
    t=0
    while True:
        t+=1
        G_plus1=update_pop(G,mutation_rate,growth_factor=growth_factor)
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
            
            