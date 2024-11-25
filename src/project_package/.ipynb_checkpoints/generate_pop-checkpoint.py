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
    
    fitnesses=np.random.rand(i) # # an array of i random numbers generated between 0 and 1
    abundances=np.random.multinomial(n,[1/i]*i) # an array of i random numbers that are roughly equal and sum to n

    G = nx.complete_bipartite_graph(1,n) # create empty graph of n strains
    for i in range(len(G.nodes)): # for each strain, add fitness and abundance data
        G.add_node(i, abundance=abundances[i],fitness=fitnesses[i])
    for start,end in G.edges:
        G.add_edge(start,end,distance=999)  
        
    return(G)

def generate_random_initial_pop_unlinked(n=100,i=20):
    ''' generate a random initial population of n individuals and i strains'''

    fitnesses=np.random.rand(i) # # an array of i random numbers generated between 0 and 1
    abundances=np.random.multinomial(n,[1/i]*i) # an array of i random numbers that are roughly equal and sum to n

    G = nx.Graph()
    G.add_nodes_from([i, {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
        
    return(G)

def generate_random_fisherlog_pop_unlinked_deprecated(i=10):
    ''' generate a random initial population of i strains whose distribution follows a Fisher-log distribution'''

    from collections import Counter
    
    a = .995 # adjusted to match Ridgeia symbs CRISPR alleles/reads
    s = np.random.logseries(a, i)
    abundances=[k for k, v in Counter(s).items() for _ in range(v)] # species_count: individual_counts in these species
    fitnesses=np.random.normal(loc=0.8, scale=0.2, size=i) #average fitness 0.8, std 0.2
    G = nx.Graph()
    G.add_nodes_from([i, {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
        
    return(G)

def generate_random_fisherlog_pop_unlinked(i=10):
    ''' generate a random initial population of i strains whose distribution follows a Fisher-log distribution. The strains are labelled "0.0.n". '''

    from collections import Counter
    
    a = .995 # adjusted to match Ridgeia symbs CRISPR alleles/reads
    s = np.random.logseries(a, i)
    abundances=[k for k, v in Counter(s).items() for _ in range(v)] # species_count: individual_counts in these species
    fitnesses=np.random.normal(loc=0.8, scale=0.2, size=i) #average fitness 0.8, std 0.2
    G = nx.Graph()
    G.add_nodes_from(['.'.join(['0.0',str(i)]), {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
        
    return(G)

def generate_random_fisherlog_pop_binomial_tree(i=10,a=0.995):
    ''' generate a random initial population of i strains whose distribution follows a Fisher-log distribution. The strains are labelled "0.0.n". '''

    from collections import Counter
    
    a = .995 # adjusted to match Ridgeia symbs CRISPR alleles/reads
    s = np.random.logseries(a, i)
    abundances=[k for k, v in Counter(s).items() for _ in range(v)] # species_count: individual_counts in these species
    fitnesses=np.random.normal(loc=0.8, scale=0.2, size=i) #average fitness 0.8, std 0.2
    print(sum(abundances))
    G = nx.binomial_tree(n=round(np.log(sum(abundances))/np.log(2)))
    nx.relabel_nodes(G, lambda x: ".".join(['0.0',str(x)]), copy=False)
    G.add_nodes_from(['.'.join(['0.0',str(i)]), {'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(fitnesses)))
    nodes_to_delete=[node for node,attr in G.nodes(data=True) if 'abundance' not in attr.keys()]
    G.remove_nodes_from(nodes_to_delete)
    nx.set_edge_attributes(G, 99, "distance")
        
    return(G)
    
    
class SymPop:

    '''
    Symbiont Population object
    '''

    
    import networkx as nx
    # from project_package.simplify import remove_empty_leaves_and_rescale_edges

    def __init__(self,G):
        self.G = G
        self.strains = [node for node,attr in G.nodes(data=True) if attr['abundance']>0]
        self.richness = len(self.strains)
        self.maxdist = sum([attr['distance'] for a,b,attr in nx.maximum_spanning_tree(G,weight='distance').edges(data=True)])
        self.mean_fitness = np.mean([attr['fitness'] for node,attr in G.nodes(data=True) if attr['abundance']>0])
        # self.heterogeneity = np.concatenate(self.get_distance_matrix()).sum()
        self.pop_size = sum([attr['abundance'] for node,attr in G.nodes(data=True)])
        # self.pi=(self.pop_size/(self.pop_size-1))*np.concatenate(self.get_pi()).sum() 
        self.gen_deg = nx.generalized_degree(G)
        self.lineages = nx.number_connected_components(G)
        self.max_degree_centrality = max([v for v in nx.degree_centrality(G).values()])
        self.Fst = 1-sum([(attr['abundance']/self.pop_size)**2 for node,attr in G.nodes(data=True) if attr['abundance']>0])
        self.Hd = self.pop_size/(self.pop_size-1) * self.Fst

    def to_ete3(self):
        import ete3

        root = min(self.G.nodes())
        subtrees = {node:ete3.Tree(name=node) for node in self.G.nodes()}
        [*map(lambda edge: subtrees[edge[0]].add_child(subtrees[edge[1]],name=edge[1],dist=str(edge[2]['distance'])), self.G.edges(data=True))]
        tree = subtrees[root]
        return(tree)

    def get_distance_matrix(self):
        header=self.strains
        p = dict(nx.shortest_path(self.G))
        print(header)
        print(p.keys())
        print(np.matrix([p[i] for i in header]))
        print(np.matrix([[p[i][j] for i in header] for j in header]))
        mat=np.matrix([[nx.path_weight(self.G, p[i][j], 'distance') for i in header] for j in header])
        return(mat)

    def get_pi(self):
        header=self.strains
        pop_size=self.pop_size
        a=np.matrix([[(self.G.nodes[i]['abundance']/pop_size)*(self.G.nodes[j]['abundance']/pop_size) for i in header] for j in header])
        b=self.get_distance_matrix()/(3500000) # estimated symbiont genome size
        return(a*b)

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


