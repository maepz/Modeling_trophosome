#### simplify the population network ####
from itertools import combinations
import networkx as nx
import numpy as np

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