#### To visualise the population network ####
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from project_package.simplify import remove_empty_leaves_and_rescale_edges

def visualize_pop_v1(G,i='0',show_empty=False):

    """
    Creates a plot of a population haplotype network. 
    OPTIONS:
    i='0' generation index to print on the graph
    show_empty=False Hide or show haplotypes with abundance=0
    """

    if show_empty==False:
        empty_leaves=[node_id  for node_id, attr in G.nodes(data=True) if (attr['abundance']==0) and (len(G.edges([node_id]))==1)]
        G=G.subgraph([n for n in G if n not in set(empty_leaves)])

    cmap = plt.get_cmap('Spectral')
    norm = mcolors.Normalize(vmin=0, vmax=1)
    node_colors = [cmap(norm(G.nodes[node]['fitness'])) for node in G.nodes()]
    
    pos = nx.spring_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=[10*attr['abundance']+1 for node_id, attr in G.nodes(data=True)])
    nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif",
                           labels={node:'\n'.join(map(str,[node,attr['abundance']])) for node,attr in G.nodes(data=True) if attr['abundance']>0})

    nx.draw_networkx_edges(G, pos, edge_color="black")

    nx.draw_networkx_edge_labels(
        G, pos, edge_labels={(start, end): attr["distance"] for start, end, attr in G.edges(data=True)}
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm,label='fitness',ax=plt.gca())
    plt.title("Population at t = "+ str(i))

def rescale(abundances):
    import numpy as np
    minsize=10 # [10^1] size of smallest node
    maxsize=10000 # [10^4] size of largest node
    pseudocount=0.01
    abundances=np.array(abundances)+pseudocount
    a=np.log10(abundances)/max(np.log10(abundances))
    b=a*(np.log10(maxsize)-np.log10(minsize))
    return(np.round(np.power(10,b))+np.log10(minsize))

def visualize_pop_v2(G,i='0',show_empty=False,view_node_labels=True,view_edge_labels=True,abundance_threshold=1):
    
    """
    Creates a plot of a population haplotype network. 
    OPTIONS:
    i='0' generation index to print on the graph
    show_empty=False Show or hide haplotypes with abundance=0; default=hide
    view_node_labels=True Show or hide node labels; default=show
    view_edge_labels=True Show or hide edges labels; default=show
    abundance_threshold=1 Abundance threshold to keep haplotype in graph. Can be absolute or relative abundance (n); default=1
    """

    from project_package.plot import rescale
    H=G.copy()
    
    # remove nodes with abundance less than x%
    if str(abundance_threshold)[-1]=='%':
        thr=np.round(sum([attr['abundance'] for node_id, attr in H.nodes(data=True)])*int(abundance_threshold[:-1])/100)
    else:
        thr=abundance_threshold
    
    H.update(nodes=[[node_id,{'abundance':0}] for node_id, attr in H.nodes(data=True) if attr['abundance'] < thr])
    
    if show_empty==False:
        H=remove_empty_leaves_and_rescale_edges(H)
              
    abundances=[attr['abundance'] for node_id, attr in H.nodes(data=True)]
    rescaled_abundances=rescale(abundances)

    # define colors
    cmap = plt.get_cmap('Spectral')
    norm = mcolors.Normalize(vmin=0, vmax=1)
    node_colors = [cmap(norm(H.nodes[node]['fitness'])) for node in H.nodes()]
    
    # draw network
    pos = nx.spring_layout(H)
    nx.draw_networkx_nodes(H, pos, node_color=node_colors, node_size=rescaled_abundances) 
    nx.draw_networkx_edges(H, pos, edge_color="black")
    
    if view_node_labels==True:
        nx.draw_networkx_labels(H, pos, font_size=12, font_family="sans-serif",
                           labels={node:'\n'.join(map(str,[node,G.nodes[node]['abundance']])) for node,attr 
                                   in H.nodes(data=True)})
    if view_edge_labels==True:        
        nx.draw_networkx_edge_labels(H, pos,edge_labels={(start, end): attr["distance"] for start, end, attr
                                                         in H.edges(data=True)})
        
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm,label='fitness',ax=plt.gca())
    plt.title("Population at t = "+ str(i))


def visualize_pop(G,i='0',show_empty=False,view_node_labels=True,view_edge_labels=True,abundance_threshold=1,transparency_threshold='100%',replace_title=False):
    
    """
    Creates a plot of a population haplotype network. 
    OPTIONS:
    i='0' generation index to print on the graph
    show_empty=False Show or hide haplotypes with abundance=0; default=hide
    view_node_labels=True Show or hide node labels; default=show
    view_edge_labels=True Show or hide edges labels; default=show
    abundance_threshold=1 Minimun abundance threshold to keep haplotype in graph. Can be absolute or relative abundance (n); default=1
    transparency_threshold=1 Maximum threshold to fill node in graph with color. Can be absolute or relative abundance (n); default='100%'
    """

    from project_package.plot import rescale
    H=G.copy()
    
    # remove nodes with abundance less than x%
    if str(abundance_threshold)[-1]=='%':
        thr=np.round(sum([attr['abundance'] for node_id, attr in H.nodes(data=True)])*int(abundance_threshold[:-1])/100)
    else:
        thr=abundance_threshold

    if str(transparency_threshold)[-1]=='%':
        fill_thr=np.round(sum([attr['abundance'] for node_id, attr in H.nodes(data=True)])*int(transparency_threshold[:-1])/100)
    else:
        fill_thr=transparency_threshold
    
    H.update(nodes=[[node_id,{'abundance':0}] for node_id, attr in H.nodes(data=True) if attr['abundance'] < thr])
    
    if show_empty==False:
        H=remove_empty_leaves_and_rescale_edges(H)

    abundances=[attr['abundance'] for node_id, attr in H.nodes(data=True)]
    rescaled_abundances=rescale(abundances)
    
    # define colors
    cmap = plt.get_cmap('Spectral')
    norm = mcolors.Normalize(vmin=0, vmax=1.5)
    node_fill = [cmap(norm(attr['fitness'])) if attr['abundance'] <= fill_thr else 'none' for node,attr in H.nodes(data=True)]
    node_outline = ['none' if attr['abundance'] <= fill_thr else 'black' for node,attr in H.nodes(data=True)]
   
    # draw network
    pos = nx.spring_layout(H)
    nx.draw_networkx_nodes(H, pos, node_color=node_fill, node_size=rescaled_abundances, edgecolors=node_outline) 
    nx.draw_networkx_edges(H, pos, edge_color="black")
    
    if view_node_labels==True:
        nx.draw_networkx_labels(H, pos, font_size=12, font_family="sans-serif",
                           labels={node:'\n'.join(map(str,[node,G.nodes[node]['abundance']])) for node,attr 
                                   in H.nodes(data=True)})
    if view_edge_labels==True:        
        nx.draw_networkx_edge_labels(H, pos,edge_labels={(start, end): attr["distance"] for start, end, attr
                                                         in H.edges(data=True)})
        
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm,label='fitness',ax=plt.gca())

    if replace_title != False:
        plt.title(replace_title)
    else:
        plt.title("Population at t = "+ str(i))