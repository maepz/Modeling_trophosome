#### simplify the population network ####
from itertools import combinations
import networkx as nx
import numpy as np

def subsample_pop(G,n):
    pop_attr=np.array([[node,attr['abundance'],attr['fitness']]for node,attr in G.nodes(data=True) if attr['abundance']>0])
    alleles = pop_attr[:,0]
    abundances = np.array(list(map(int,pop_attr[:,1])))
    fitnesses = np.array(list(map(float,pop_attr[:,2])))

    tot=sum(abundances)
    weights=np.multiply(abundances,1/tot)
    new_pop_abundances = np.random.multinomial(n, weights)
    new_pop_alleles = [alleles[i] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    new_G=nx.Graph(G.subgraph(new_pop_alleles))
    
    adj=[[alleles[i],{'abundance':new_pop_abundances[i],'fitness':fitnesses[i]}] for i in range(len(new_pop_abundances)) if new_pop_abundances[i] > 0]
    new_G.update(nodes=adj)

    return(new_G)
    
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


def remove_empty_nodes_and_reconnect(G):
    ''' DEPRECATED. This version is not good because it loops through the nodes. Instead, the long branches should be extracted first.'''
    nodes_to_remove = [node for node,attr in G.nodes(data=True) if attr['abundance']==0]
    newG=G.copy()
    
    for node in nodes_to_remove:
        old_edges=G.edges([node])
        new_edges=list(combinations([other_node for node, other_node in old_edges],2))
        new_edges_with_weights=[(edge[0],edge[1],{'distance':nx.shortest_path_length(G, source=edge[0], target=edge[1], weight='distance')}) for edge in new_edges]
        newG.add_edges_from(new_edges_with_weights)
        newG.remove_edges_from(old_edges)
    newG.remove_nodes_from(nodes_to_remove)
    
    return(newG)


def remove_empty_leaves_and_rescale_edges(G):

    '''
    This function simplify the haplotype network by removing haplotypes that have null abundances and rescaling the distance between the remaining haplotypes
    '''
    newG=G.copy()
    
    ## remove empty leaves
    
    nodes_to_remove = [node for node,attr in G.nodes(data=True) if attr['abundance']==0 and G.degree(node)<2]
    edges_to_remove = G.edges(nodes_to_remove)
    
    newG.remove_nodes_from(nodes_to_remove)
    newG.remove_edges_from(edges_to_remove)
    
    ## find and collapse chain of edges between empty nodes
    
    nodes_to_remove = [node for node,attr in newG.nodes(data=True) if attr['abundance']==0 and newG.degree(node)==2]
    edges_to_remove = []
    new_edges_with_weights =[]
    
    subgraphs=[newG.subgraph(c).copy() for c in nx.connected_components(newG.subgraph(nodes_to_remove))]

    for sg in subgraphs:
        old_edges=newG.edges(sg.nodes())
        edges_to_remove+=old_edges
        
        new_edge=[node for node in np.array(list(old_edges)).flatten() if node not in sg.nodes()]
        new_edges_with_weights+=[[new_edge[0],new_edge[1],{'distance':nx.shortest_path_length(newG, source=new_edge[0], target=new_edge[1], weight='distance')}]]

    newG.remove_nodes_from(nodes_to_remove)
    newG.remove_edges_from(edges_to_remove)
    newG.add_edges_from(new_edges_with_weights)
    
    return(newG)

def remove_empty_leaves(G):

    '''
    This function simplify the haplotype network by removing leafs haplotypes that have null abundances.
    '''
    newG=G.copy()
    
    ## remove empty leaves
    
    nodes_to_remove = [node for node,attr in G.nodes(data=True) if attr['abundance']==0 and G.degree(node)<2]
    edges_to_remove = G.edges(nodes_to_remove)
    
    newG.remove_nodes_from(nodes_to_remove)
    newG.remove_edges_from(edges_to_remove)
        
    return(newG)

def merge_graphs(Graph_list):
    
    '''This function takes multiple graphs (from a list of networkx objects) and merges them into a single one. Combining common nodes and edges, summing the Node abundances across graphs.'''
    
    merged_Graph=nx.compose_all(Graph_list)

    ### sum abundances from shared nodes and update merged graph
    dic_list=[nx.get_node_attributes(Graph_list[i],'abundance') for i in range(len(Graph_list))]
    from functools import reduce
    def reducer(accumulator, element):
        for key, value in element.items():
            accumulator[key] = accumulator.get(key, 0) + value
        return accumulator
    node_data=reduce(reducer, dic_list, {})
    nx.set_node_attributes(merged_Graph, node_data, "abundance")

    return(merged_Graph)
