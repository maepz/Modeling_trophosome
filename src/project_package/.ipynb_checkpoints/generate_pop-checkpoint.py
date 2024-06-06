#### To generate the initial population ####

import networkx as nx
import numpy as np

def generate_initial_pop(fitnesses,abundances):
    ''' generate an initial population of len(fitnesses) strains with fitness and abundance values provided in lists'''
    
    # fitnesses=np.insert(fitnesses,0, 0)
    fitnesses=[0]+fitnesses
    abundances=[0]+abundances
    # abundances=np.insert(abundances,0, 0)

    G = nx.complete_bipartite_graph(1,len(fitnesses)-1) # create empty graph of n strains
    for i in range(len(G.nodes)): # for each strain, add fitness and abundance data
        G.add_node(i, abundance=abundances[i],fitness=fitnesses[i])
    for start,end in G.edges:
        G.add_edge(start,end,distance=999)  
        
    return(G)


def generate_initial_pop_unlinked(fitnesses,abundances):
    ''' generate an initial population of len(fitnesses) strains with fitness and abundance values provided in lists'''
    
    G = nx.Graph()
    G.add_nodes_from([i, {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
            
    return(G)

def generate_random_initial_pop(n):
    ''' generate a random initial population of n strains'''
    
    fitnesses=np.insert(np.random.rand(n), 0, 0) # # an array of n random numbers generated between 0 and 1
    abundances=np.insert(np.ceil(np.random.rand(n)*100).astype(int),0,0) # an array of n random numbers generated between 1 and 100

    # print(fitnesses)
    # print(abundances)

    G = nx.complete_bipartite_graph(1,n) # create empty graph of n strains
    for i in range(len(G.nodes)): # for each strain, add fitness and abundance data
        G.add_node(i, abundance=abundances[i],fitness=fitnesses[i])
    for start,end in G.edges:
        G.add_edge(start,end,distance=999)  
        
    return(G)

def generate_random_initial_pop_unlinked(n=100,i=20):
    ''' generate a random initial population of n individuals and i strains'''

    fitnesses=np.random.rand(n) # # an array of n random numbers generated between 0 and 1
    abundances=np.ceil(np.random.rand(n)*100).astype(int) # an array of n random numbers generated between 1 and 100

    G = nx.Graph()
    G.add_nodes_from([i, {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
        
    return(G)

def generate_random_fisherlog_pop_unlinked(i=10):
    ''' generate a random initial population of i strains whose distribution follows a Fisher-log distribution'''

    from collections import Counter
    
    a = .995 # adjusted to match Ridgeia symbs CRISPR alleles/reads
    s = np.random.logseries(a, i)
    abundances=[k for k, v in Counter(s).items() for _ in range(v)] # species_count: individual_counts in these species
    fitnesses=np.random.normal(loc=0.8, scale=0.2, size=i) #average fitness 0.8, std 0.2
    G = nx.Graph()
    G.add_nodes_from([i, {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
        
    return(G)

class SymPop:

    '''
    Symbiont Population object
    '''

    
    import networkx as nx
    # from project_package.simplify import remove_empty_leaves_and_rescale_edges

    def __init__(self, G):
        self.G = G
        self.richness = len([node for node,attr in G.nodes(data=True) if attr['abundance']>0])
        self.maxdist = sum([attr['distance'] for a,b,attr in nx.maximum_spanning_tree(G,weight='distance').edges(data=True)])
        self.mean_fitness = np.mean([attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0])
        self.heterogeneity = sum([attr['distance'] for a,b, attr in G.edges(data=True)])
        self.pop_size = sum([attr['abundance'] for node,attr in G.nodes(data=True)])
        self.gen_deg = nx.generalized_degree(G)
        self.lineages = nx.number_connected_components(G)
        self.max_degree_centrality = max([v for v in nx.degree_centrality(G).values()])
        self.Fst = 1-sum([(attr['abundance']/self.pop_size)**2 for node,attr in G.nodes(data=True) if attr['abundance']>0])

    def to_ete3(self):
        import ete3

        root = min(self.G.nodes())
        subtrees = {node:ete3.Tree(name=node) for node in self.G.nodes()}
        [*map(lambda edge: subtrees[edge[0]].add_child(subtrees[edge[1]],name=edge[1],dist=str(edge[2]['distance'])), self.G.edges(data=True))]
        tree = subtrees[root]
        return(tree)

    def visualize_pop(self,i='0',show_empty=False,view_node_labels=True,view_edge_labels=True,abundance_threshold=1):
        from project_package.plot import visualize_pop
        show_empty=show_empty
        view_node_labels=view_node_labels
        view_edge_labels=view_edge_labels
        abundance_threshold=abundance_threshold
        visualize_pop(self.G,show_empty=show_empty,
                      view_node_labels=view_node_labels,view_edge_labels=view_edge_labels,
                      abundance_threshold=abundance_threshold)

def remove_empty_leaves_and_rescale_edges(self):
        from project_package.simplify import remove_empty_leaves_and_rescale_edges
        return(remove_empty_leaves_and_rescale_edges(self))