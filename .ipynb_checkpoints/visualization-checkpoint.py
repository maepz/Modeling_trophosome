#### To visualise the population network ####

def visualize_pop(G,i='0'):
    
    cmap = plt.get_cmap('Spectral')
    norm = mcolors.Normalize(vmin=0, vmax=1)
    node_colors = [cmap(norm(G.nodes[node]['fitness'])) for node in G.nodes()]

    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=[10*attr['abundance'] for node_id, attr in G.nodes(data=True)])
    nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif",
                           labels={node:'\n'.join(map(str,[node,attr['abundance']])) for node,attr in G.nodes(data=True) if attr['abundance']>0}
                           )
    # if len(G.edges)>0:
    nx.draw_networkx_edges(G, pos, edge_color="black")

    nx.draw_networkx_edge_labels(
        G, pos, edge_labels={(start, end): attr["distance"] for start, end, attr in G.edges(data=True)}
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm,label='fitness',ax=plt.gca())
    plt.title("Population at t = "+ str(i))