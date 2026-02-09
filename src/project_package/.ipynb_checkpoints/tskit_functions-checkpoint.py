import tskit
import numpy as np
import json
import msprime
import pandas as pd
import math
import networkx as nx
from collections import defaultdict, deque

# def simulate_ancestry_with_mutations(new_avail_id,n, baseline=1, random_seed=None, sequence_length=1, mu=1e-9):
#     ts = msprime.sim_ancestry(n, random_seed=seed, sequence_length=1)
#     mts = msprime.sim_mutations(
#         ts, 
#         rate=mu, 
#         model=msprime.SLiMMutationModel(type=1),
#         random_seed=seed)

#     if len([var for var in mts.variants()])==0:
#         print('no variants')
#         return([],[])

#     final_intrahost_strains=+[var.counts() for var in mts.variants()][0]
#     final_intrahost_strains_fitness_dic=fitnesses_by_genotype(mts, )
#     print('final_intrahost_strains.keys()',final_intrahost_strains.keys())

#     parent_strain_final_count=0
#     try:
#         parent_strain_final_count=final_intrahost_strains.pop('')
#     except KeyError:
#         print(list(final_intrahost_strains_fitness_dic.keys()))
#         ml = {m.id: m.derived_state for m in mts.mutations()}
#         plt.figure()
#         mts.draw_svg(mutation_labels=ml)
#         return([],[])
        
#     alleles=np.array(list(final_intrahost_strains.keys()))
#     abundances=list(final_intrahost_strains.values())
#     distances=[sum(tok.isdigit() for tok in s.split(',')) for s in final_intrahost_strains.keys()]
#     fitnesses=[final_intrahost_strains_fitness_dic[allele] for allele in alleles]
        
#     # Update pop with new alleles    
#     new_avail_prefix='.'.join(new_avail_id.split('.')[0:-1])
#     new_avail_idx=int(new_avail_id.split('.')[-1])
#     new_alleles_ids=['.'.join([new_avail_prefix,str(i)]) for i in range(new_avail_idx,new_avail_idx+len(final_intrahost_strains))]
#     new_avail_id='.'.join([new_avail_prefix,str(new_avail_idx+len(final_intrahost_strains))])

                
#     # new_alleles_parent_indices=np.repeat(mutated_alleles, mutated_ind[mutated_alleles]).astype(int)
#     # mutated_alleles_parent_fitness=fitnesses[new_alleles_parent_indices]
#     # new_alleles_fitnesses=np.sum([mutated_alleles_parent_fitness,np.random.normal(-0.01, 0.01,size=mutated_ind_count)], axis=0) # Assign fitness to the new alleles. The fitness of the new alleles is that of the parent allele +/- a selection coefficient sampled from a normal distribution with mean -0.01 and std 0.01
    
#     new_alleles_nodes=[[new_alleles_ids[i],{'abundance':abundances[i],'fitness':fitnesses[i]}] for i in range(len(alleles))]
#     # new_alleles_edges=np.array([alleles[new_alleles_parent_indices],new_alleles_ids,[{"distance":1}] * len(new_alleles_ids)]).transpose()
#     new_alleles_edges=[]
#     return(new_alleles_nodes,new_alleles_edges)

# def relabel_networkx_graph(G,label_dic):
    


def tstree_to_networkx_graph(ts,Trim=True):
    edges=ts.tables.edges
    nodes=ts.tables.nodes

    G = nx.from_pandas_edgelist(
        pd.DataFrame({'source': edges.parent, 'target': edges.child, 'left': edges.left, 'right': edges.right,'distance':[d.metadata["distance"] for d in edges]}),
        edge_attr=True,
        create_using=nx.Graph
    )

    # nx.set_edge_attributes(G, {i: {'distance':edges[i].metadata["distance"]} for i in range(len(edges))})
    nx.set_node_attributes(G, {i: {'abundance':nodes[i].metadata["abundance"], 'fitness': nodes[i].metadata["fitness"]} for i in range(len(nodes))})



    if Trim==True:
        """
        Merge nodes connected by edges with distance == 0.
        Representative choice rule:
          - If any node in the component has abundance > 0, choose the one with the
            largest abundance; ties broken by larger node id.
          - If all nodes have abundance == 0, choose the node with the larger id.
        After relabeling, remove all edges with distance == 0 and any self-loops.
        Node attributes are exactly those of the representative (no 'contraction' info).
        """
        #### V1
        # newG=G.copy()
        # edges_to_remove = [(node_1,node_2) for node_1,node_2,attr in G.edges(data=True) if attr['distance']==0]
        # newG.remove_edges_from(edges_to_remove)

        # return(newG)
        
        ##### V2
        # edges_to_remove = [(node_1,node_2) for node_1,node_2,attr in G.edges(data=True) if attr['distance']==0]
        # while len(edges_to_remove)>0:
        #     edge=edges_to_remove.pop()
        #     node1=sorted(edge)[0]
        #     node2=sorted(edge)[1]
        #     G=nx.contracted_nodes(G, node1, node2,self_loops=False)

        #### V3
        # Build an undirected helper graph of zero-distance edges
        U = nx.Graph()
        U.add_nodes_from(G.nodes())
        if G.is_multigraph():
            for u, v, k, data in G.edges(keys=True, data=True):
                if data.get("distance") == 0:
                    U.add_edge(u, v)
        else:
            for u, v, data in G.edges(data=True):
                if data.get("distance") == 0:
                    U.add_edge(u, v)
    
        # Choose representatives per component using your rule
        def abund(n):
            return G.nodes[n].get("abundance", 0)
    
        mapping = {}
        for comp in nx.connected_components(U):
            comp = list(comp)
            # Prefer nodes with abundance > 0, pick largest abundance then larger id
            positives = [n for n in comp if abund(n) > 0]
            if positives:
                rep = max(positives, key=lambda n: (abund(n), n))
            else:
                # All zero: keep the node with the larger id
                rep = max(comp)
            for n in comp:
                mapping[n] = rep
    
        # Relabel nodes to their representatives (merges nodes); then clean up edges
        H = nx.relabel_nodes(G, mapping, copy=True)
    
        # Ensure node attributes come from the representative exactly
        reps = set(mapping.values())
        for rep in reps:
            if rep in G:
                # overwrite any merged attributes with those of the representative
                H.nodes[rep].clear()
                H.nodes[rep].update(G.nodes[rep])
    
        # Remove zero-distance edges and self-loops
        if H.is_multigraph():
            H.remove_edges_from([(u, v, k) for u, v, k, d in H.edges(keys=True, data=True)
                                 if d.get("distance") == 0])
            H.remove_edges_from(list(H.selfloop_edges(keys=True)))
        else:
            H.remove_edges_from([(u, v) for u, v, d in H.edges(data=True)
                                 if d.get("distance") == 0])
            H.remove_edges_from(list(nx.selfloop_edges(H)))
    
        return( H)

    return(G)


def cell_coalescent_tree_with_mutations(n, mu, fitness, tree_seed, mutation_seed, fitness_seed):

    ts = msprime.sim_ancestry(n, random_seed=tree_seed, sequence_length=1)
    mts = msprime.sim_mutations(
        ts, 
        rate=mu, 
        model=msprime.SLiMMutationModel(type=1),
        random_seed=mutation_seed)

    mts_unique=prune_unique_genotypes_add_metadata(mts, n, s_mean=-0.1, s_std=0.1, baseline=fitness, seed=fitness_seed)
    
    return(mts_unique)

def prune_unique_genotypes_add_metadata(ts, n, s_mean=-0.1, s_std=0.1, baseline=1.0, seed=None):
    """
    Returns a new TreeSequence pruned to one sample per unique genotype and annotated with:
    - mutation metadata: {"s": float} with s ~ Normal(s_mean, s_std)
    - node metadata: {"fitness": float, "abundance": int, "mutation_count": int}
    - edge metadata: {"distance": int} = number of mutations on that branch
    Assumes haploid samples. Handles SLiMMutationModel stacking at each site.
    """
    rng = np.random.default_rng(seed)
    
    # 0) Define schema

    import tskit

    tables = ts.dump_tables()

    tables.edges.metadata_schema = tskit.MetadataSchema({
    "codec": "json",
    "schema": {
    "type": "object",
    "properties": {
    "distance": {"type": "integer"},
    },
    "additionalProperties": False
    }
    })

    tables.nodes.metadata_schema = tskit.MetadataSchema({
    "codec": "json",
    "schema": {
    "type": "object",
    "properties": {
    "fitness": {"type": "float"},
    "abundance":{"type": "integer"},
    "mutation_count": {"type": "integer"}
    },
    "additionalProperties": False
    }
    })

    tables.mutations.metadata_schema = tskit.MetadataSchema({
    "codec": "json",
    "schema": {
    "type": "object",
    "properties": {
    "s": {"type": "float"}
    },
    "additionalProperties": False
    }
    })

    # 1) Annotate mutations with {"s": ...}
    muts = tables.mutations
    # print('muts',muts)
    n_muts = muts.num_rows
    # print('n_muts',n_muts)
    # print(tables.nodes)

    if n_muts ==0: ## If there is no mutations, return a simple tree

        tables.nodes.clear()
        tables.edges.clear()
        tables.sites.clear()
        tables.mutations.clear()

        tables.nodes.metadata_schema = tskit.MetadataSchema({
            "codec": "json",
            "schema": {
            "type": "object",
            "properties": {
            "fitness": {"type": "float"},
            "abundance":{"type": "integer"},
            "mutation_count": {"type": "integer"}
            },
            "additionalProperties": False
            }
            })

        n0 = ts.node(0)
        # Add a root with an older time
        root_time = float(n0.time) + 1.0
        root_id = tables.nodes.add_row(
            time=root_time,
            flags=0,
            population=tskit.NULL,
            individual=tskit.NULL,
            metadata={"fitness": float(baseline), "abundance": 0, "mutation_count": 0},  # no node metadata required for the root
        )
        
        # Add the single sample first so it is id 0
        sample_id = tables.nodes.add_row(
            time=n0.time,
            flags=tskit.NODE_IS_SAMPLE,  # ensure it’s a sample
            population=n0.population,
            individual=n0.individual,
            metadata={"fitness": float(baseline), "abundance": int(n), "mutation_count": 0})
        
        
        
        # Connect root -> sample across the whole genome; set distance=0 to satisfy edge schema
        tables.edges.add_row(
            left=0.0,
            right=ts.sequence_length,
            parent=root_id,
            child=sample_id,
            metadata={"distance": 0},
        )
        
        tables.sort()
        return tables.tree_sequence()

        
        tables.nodes.add_row(
            time=ts1.node(0).time,
            flags=ts1.node(0).flags | tskit.NODE_IS_SAMPLE,  # ensure it’s a sample
            population=ts1.node(0).population,
            individual=ts1.node(0).individual,
            metadata={ "fitness": baseline, "abundance": n, "mutation_count": 0})
        # tables.sort()
        # return tables.tree_sequence()
        
    else:
        s_vals = rng.normal(loc=s_mean, scale=s_std, size=n_muts)
        original_muts = list(ts.mutations())
        muts.clear()
        # tskit 1.0.0 does not have mutation time; guard in case
        has_time = hasattr(original_muts[0], "time")
        for i, m in enumerate(original_muts):
            if has_time:
                muts.add_row(site=m.site, node=m.node, derived_state=m.derived_state,
                             parent=m.parent, time=m.time, metadata={"s":float(s_vals[i])})
            else:
                muts.add_row(site=m.site, node=m.node, derived_state=m.derived_state,
                             parent=m.parent, metadata={"s":float(s_vals[i])})
    ts_s = tables.tree_sequence()
    
    # 2) Group haploid samples by unique multilocus genotype
    samples = ts_s.samples()
    num_samples = len(samples)
    # Stream genotypes per site to avoid full matrix
    genotype_components = [[] for _ in range(num_samples)]
    for var in ts_s.variants():
        g = var.genotypes
        for i in range(num_samples):
            genotype_components[i].append(int(g[i]))
    genotype_groups = {}
    for i in range(num_samples):
        key = tuple(genotype_components[i])  # empty tuple if no sites
        genotype_groups.setdefault(key, []).append(i)
    
    # Pick one representative per group; abundance = group size
    rep_indices = []
    abundance_by_rep_node = {}
    for idxs in genotype_groups.values():
        rep_i = min(idxs, key=lambda j: samples[j])  # stable representative choice
        rep_indices.append(rep_i)
        abundance_by_rep_node[samples[rep_i]] = len(idxs)

    rep_nodes = [samples[i] for i in rep_indices]
    
    # 3) Pre-annotate abundance for representative nodes in tables, then simplify
    nt = tables.nodes
    num_nodes = nt.num_rows
    flags = nt.flags
    time = nt.time
    population = nt.population
    individual = getattr(nt, "individual", np.full(num_nodes, tskit.NULL, dtype=np.int32))
    # Rebuild node table with abundance metadata
    nt.clear()
    for nid in range(num_nodes):
        nt.add_row(flags=flags[nid], time=time[nid], population=population[nid],
                   individual=individual[nid], metadata={"abundance": int(abundance_by_rep_node.get(nid, 0))})

    
    ts_pruned = tables.tree_sequence().simplify(samples=rep_nodes, keep_unary=False, filter_sites=False)
    
    # 4) Compute per-node fitness and cumulative mutation count on the pruned TS
    node_s_accum = np.zeros(ts_pruned.num_nodes, dtype=float)
    node_mut_accum = np.zeros(ts_pruned.num_nodes, dtype=np.int32)

    for site in ts_pruned.sites():
        site_muts = list(site.mutations)
        if not site_muts:
            continue
    
    try:
        tree = ts_pruned.at(site.position)
    except:
        print(tables)
        print(ts_pruned)
        print(ts_pruned.tables.edges)
        print(ts_pruned.tables.nodes)
        print(ts_pruned.tables.sites)

    roots = tree.roots if hasattr(tree, "roots") else ([tree.root] if tree.root != tskit.NULL else [])
    root_set = set(roots)
    
    # Per-site effects on child nodes
    edge_s = np.zeros(ts_pruned.num_nodes, dtype=float)
    edge_cnt = np.zeros(ts_pruned.num_nodes, dtype=np.int32)
    
    for m in site_muts:
        # Skip mutations placed on roots so root stays at baseline fitness and count starts at 0
        if m.node not in root_set:
            edge_s[m.node] += float(m.metadata["s"])
            edge_cnt[m.node] += 1
    
    # Propagate s and counts from each root to tips
    for r in roots:
        stack = [(r, 0.0, 0)]
        while stack:
            u, s_contrib, c_contrib = stack.pop()
            node_s_accum[u] += s_contrib
            node_mut_accum[u] += c_contrib
            for v in tree.children(u):
                stack.append((v, s_contrib + edge_s[v], c_contrib + edge_cnt[v]))
    node_fitness = baseline + node_s_accum
    
    # 5) Final node metadata: {"fitness", "abundance", "mutation_count"} using cumulative counts
    tables2 = ts_pruned.dump_tables()
    nt2 = tables2.nodes
    num_nodes2 = nt2.num_rows
    flags2 = nt2.flags
    time2 = nt2.time
    population2 = nt2.population
    individual2 = getattr(nt2, "individual", np.full(num_nodes2, tskit.NULL, dtype=np.int32))
    abundance = dict(zip(range(nt2.num_rows),[row.metadata["abundance"] for row in nt2]))
    mutation_count = node_mut_accum.astype(np.int32)
    
    # Rebuild node table with final metadata
    nt2.clear()
    for nid in range(num_nodes2):
        nt2.add_row(
        flags=flags2[nid],
        time=time2[nid],
        population=population2[nid],
        individual=individual2[nid],
        metadata={
        "fitness": float(node_fitness[nid]),
        "abundance": int(abundance[nid]),
        "mutation_count": int(mutation_count[nid]),
        },
        )    
    
    # 6) Edge metadata: {"distance": number of mutations on that branch within [left, right)}
    et = tables2.edges
    num_edges = et.num_rows
    left = et.left
    right = et.right
    parent = et.parent
    child = et.child
    
    edge_meta = [""] * num_edges
    if num_edges > 0 and tables2.mutations.num_rows > 0:
        m_nodes = tables2.mutations.node
        m_sites = tables2.mutations.site
        site_pos = tables2.sites.position
        m_pos = site_pos[m_sites]
    
        order = np.lexsort((m_pos, m_nodes))
        nodes_sorted = m_nodes[order]
        pos_sorted = m_pos[order]
    
        unique_nodes, starts, counts = np.unique(nodes_sorted, return_index=True, return_counts=True)
        nrange_start = np.full(tables2.nodes.num_rows, -1, dtype=np.int64)
        nrange_count = np.zeros(tables2.nodes.num_rows, dtype=np.int64)
        nrange_start[unique_nodes] = starts
        nrange_count[unique_nodes] = counts
    
        for i in range(num_edges):
            c = child[i]
            s = nrange_start[c]
            cnt = 0
            if s >= 0:
                seg = pos_sorted[s : s + nrange_count[c]]
                lo = np.searchsorted(seg, left[i], side="left")
                hi = np.searchsorted(seg, right[i], side="left")
                cnt = int(hi - lo)
            edge_meta[i] = {"distance": cnt} 
    else:
        for i in range(num_edges):
            edge_meta[i] = {"distance": 0}
    
    # Rebuild edge table with metadata
    orig_left = np.array(left, copy=True)
    orig_right = np.array(right, copy=True)
    orig_parent = np.array(parent, copy=True)
    orig_child = np.array(child, copy=True)
    
    et.clear()
    for i in range(num_edges):
        et.add_row(left=orig_left[i], right=orig_right[i], parent=orig_parent[i],
                   child=orig_child[i], metadata=edge_meta[i])
    
    return(tables2.tree_sequence())


def NodeTable_to_pandasDF(nt): #nt=mts.tables.nodes
    df = pd.DataFrame({
    "id": np.arange(nt.num_rows, dtype=np.int32),
    "flags": np.array(nt.flags, copy=False),
    "time": np.array(nt.time, copy=False),
    "population": np.array(nt.population, copy=False),
    "individual": np.array(getattr(nt, "individual",
    np.full(nt.num_rows, -1, dtype=np.int32)),
    copy=False),
    })
    md_rows = [row.metadata for row in nt]
    
    if md_rows and isinstance(md_rows[0], dict):
        md_df = pd.DataFrame(md_rows).fillna(pd.NA)
        df = pd.concat([df, md_df], axis=1)
    else:
        # If bytes, try JSON decode; otherwise keep raw bytes in a single column
        try:
            decoded = [json.loads(m) if isinstance(m, (bytes, bytearray)) and len(m) > 0 else {}
            for m in md_rows]
            md_df = pd.DataFrame(decoded).fillna(pd.NA)
            df = pd.concat([df, md_df], axis=1)
        except Exception:
            df["metadata"] = md_rows
    return(df)

def merge_by_root_and_combine_ancestrals_v1(
ts_list,
tree_index=0,
prefixes=None,
sequence_length=1.0,
root_time=None,
label_mode="all",
root_metadata=None,
site_ancestral_state=None,
site_metadata=None,
preserve_branch_lengths_when_possible=True,
shift_mutation_times=True,
remove_move_leaf_sample_flag=True,
epsilon=1e-9,
):
    """
    Merge selected trees from multiple TreeSequences by the root (single SLiM site),
    then prune & regraft duplicate terminal lineages:
    - keep_leaf = first terminal sample with metadata[abundance] > 0 and metadata[mutation_count] == 0,
    - move_leafs = all subsequent terminal samples with the same condition.
    For each move_leaf:
    * Regraft all other children of its parent (excluding move_leaf) under keep_leaf's parent,
    * Remove edges (root -> p_move) and (p_move -> move_leaf),
    * Add move_leaf's abundance into keep_leaf's abundance (integer).
    Optionally clear move_leaf sample flags.
    
    Returns:
      merged_ts: TreeSequence after merge+prune/regraft
      node_labels: dict mapping output node id -> label
    """
    if len(ts_list)==1:
        return(ts_list[0])
    if prefixes is None:
        prefixes = [f"{i}_" for i in range(len(ts_list))]
    if len(prefixes) != len(ts_list):
        raise ValueError("prefixes length must match number of inputs")
    
    # Helper: position a Tree at the left coordinate of the idx-th tree
    def get_tree(ts, idx):
        if idx < 0 or idx >= ts.num_trees:
            raise IndexError(f"tree_index {idx} out of range for input with {ts.num_trees} trees")
        lefts = list(ts.breakpoints())
        return ts.at(lefts[idx])
    
    # Collect trees and root times
    trees = []
    roots = []
    root_times = []
    for i, ts in enumerate(ts_list):
        t = get_tree(ts, tree_index)
        rts = t.roots if hasattr(t, "roots") else ([t.root] if t.root != tskit.NULL else [])
        if len(rts) != 1:
            raise ValueError(f"Input {i} has {len(rts)} roots; expected 1")
        r = rts[0]
        trees.append(t)
        roots.append(r)
        root_times.append(ts.node(r).time)
    
    # Choose common root time and per-tree time shifts to preserve branch lengths
    T = max(root_times) if root_time is None else float(root_time)
    shifts = [T - rt for rt in root_times]
    
    tables = tskit.TableCollection(sequence_length=float(sequence_length))
    
    # Set your metadata schemas
    tables.edges.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {
            "type": "object",
            "properties": {"distance": {"type": "integer"}},
            "additionalProperties": False,
        },
    })
    tables.nodes.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {
            "type": "object",
            "properties": {
                "fitness": {"type": "float"},
                "abundance": {"type": "integer"},
                "mutation_count": {"type": "integer"},
            },
            "additionalProperties": False,
        },
    })
    tables.mutations.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {
            "type": "object",
            "properties": {"s": {"type": "float"}},
            "additionalProperties": False,
        },
    })
    tables.sites.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {"type": "object", "properties": {}, "additionalProperties": False},
    })
    
    # Single SLiM site first: id=0, position=0
    if site_ancestral_state is None:
        site_ancestral_state = ""
        for ts in ts_list:
            if ts.num_sites > 0:
                site_ancestral_state = ts.tables.sites[0].ancestral_state
                break
    if site_metadata is None:
        for ts in ts_list:
            if ts.num_sites > 0:
                md = ts.tables.sites[0].metadata
                if md not in (None, {}, b"", "", []):
                    site_metadata = md
                    break
    if site_metadata is None:
        site_metadata = {}
    site0_id = tables.sites.add_row(position=0.0, ancestral_state=site_ancestral_state, metadata=site_metadata)
    assert site0_id == 0, "Site id is expected to be 0"
    
    # Synthetic root node
    root_out_md = {"fitness": np.nan,
                "abundance": 0,
                "mutation_count": 0} if root_metadata is None else dict(root_metadata)
    root_out = tables.nodes.add_row(time=T, flags=0, population=tskit.NULL, individual=tskit.NULL, metadata=root_out_md)
    
    # Map input node IDs to output node IDs and build labels
    node_maps = []
    node_labels = {}
    for i, ts in enumerate(ts_list):
        m = np.full(ts.num_nodes, -1, dtype=np.int32)
        prefix = prefixes[i]
        shift = shifts[i]
        is_sample = np.zeros(ts.num_nodes, dtype=bool)
        is_sample[ts.samples()] = True
        for u in range(ts.num_nodes):
            n = ts.node(u)
            new_id = tables.nodes.add_row(
                time=n.time + shift,
                flags=n.flags,
                population=tskit.NULL,  # keep NULL; copy populations/individuals if needed
                individual=tskit.NULL,
                metadata=n.metadata,    # dict; schema handles encoding
            )
            m[u] = new_id
            if label_mode == "all" or (label_mode == "leaves" and is_sample[u]):
                node_labels[new_id] = f"{prefix}{u}"
        node_maps.append(m)
    
    # Add edges over the full interval; reparent root children to synthetic root
    SL = float(sequence_length)
    for i, (ts, tree, r_in) in enumerate(zip(ts_list, trees, roots)):
        m = node_maps[i]
        for u in tree.nodes():
            p = tree.parent(u)
            if p == tskit.NULL:
                continue
            # Try to fetch input edge metadata for (p -> u) in this tree
            e_id = getattr(tree, "edge", lambda _: tskit.NULL)(u)
            md = ts.tables.edges[e_id].metadata if e_id != tskit.NULL else None
            parent_out = root_out if p == r_in else m[p]
            child_out = m[u]
            tables.edges.add_row(left=0.0, right=SL, parent=parent_out, child=child_out, metadata=md)
    
    # Copy mutations: remap node to merged ID, site -> 0, shift time; copy metadata
    for i, ts in enumerate(ts_list):
        m_map = node_maps[i]
        has_mut_time = hasattr(ts.tables.mutations, "time") and len(ts.tables.mutations.time) == ts.num_mutations
        shift = shifts[i]
        mut_id_map = np.full(ts.num_mutations, -1, dtype=np.int32)
        pending = list(range(ts.num_mutations))
        while pending:
            progressed = False
            next_pending = []
            for old_mid in pending:
                m = ts.tables.mutations[old_mid]
                parent_ok = (m.parent == tskit.NULL) or (mut_id_map[m.parent] != -1)
                if parent_ok:
                    new_parent = tskit.NULL if m.parent == tskit.NULL else mut_id_map[m.parent]
                    new_mid = tables.mutations.add_row(
                        site=0,
                        node=m_map[m.node],
                        derived_state=m.derived_state,
                        parent=new_parent,
                        time=(m.time + shift) if has_mut_time and m.time is not None else None,
                        metadata=m.metadata,
                    )
                    mut_id_map[old_mid] = new_mid
                    progressed = True
                else:
                    next_pending.append(old_mid)
            if not progressed:
                raise ValueError("Could not resolve mutation.parent references; cycle or missing parent?")
            pending = next_pending
    
    # Sort and build a temporary merged TreeSequence to find leaves/parents
    tables.sort()
    merged_ts = tables.tree_sequence()
    
    # Identify keep_leaf and move_leafs in the merged tree
    tree = merged_ts.first()
    root = tree.root
    if root == tskit.NULL:
        raise ValueError("Merged tree has no root")
    
    def is_terminal_sample(u):
        return ((merged_ts.node(u).flags & tskit.NODE_IS_SAMPLE) != 0) and (tree.num_children(u) == 0)
    
    def md_ok(u):
        md = merged_ts.node(u).metadata or {}
        return (md.get('abundance', 0) > 0) and (md.get('mutation_count', 0) == 0)
    
    leaves_in_order = list(tree.leaves(root))
    candidates = [u for u in leaves_in_order if is_terminal_sample(u) and md_ok(u)]
    if len(candidates) <= 1:
        # Nothing to merge further; return the merged tree and labels
        return( merged_ts, node_labels)
    
    keep_leaf = candidates[0]
    move_leafs = candidates[1:]
    
    # Validate parents under root
    p_keep = tree.parent(keep_leaf)
    # print('keep_leaf',keep_leaf,'p_keep',p_keep, 'root',root)
    # print(merged_ts.tables.nodes)
    # print(merged_ts.tables.edges)
    # print('tree.parent(p_keep)',tree.parent(p_keep))
    # print(p_keep)
    if p_keep == tskit.NULL or tree.parent(p_keep) != root:
        raise ValueError("keep_leaf must be under an internal parent that is a direct child of the root")
        # print("keep_leaf must be under an internal parent that is a direct child of the root")
        # return(merged_ts)

    
    move_by_parent = {}
    p_moves = set()
    for u in move_leafs:
        p = tree.parent(u)
        if p == tskit.NULL or tree.parent(p) != root:
            raise ValueError(f"move_leaf {u} must be under an internal parent that is a child of the root")
        move_by_parent[p] = u
        p_moves.add(p)
    
    # Modify tables (dump) to regraft subtrees and remove branches
    tables = merged_ts.dump_tables()
    
    # children map from current edges
    children_of = defaultdict(set)
    for e in tables.edges:
        children_of[e.parent].add(e.child)
    
    # For each move parent, determine children to reattach under p_keep
    roots_to_move = {}
    for p in p_moves:
        mv = move_by_parent[p]
        roots_to_move[p] = [c for c in children_of[p] if c != mv]
    
    # Shift times for moved subtrees when p_move older than p_keep
    t_keep = merged_ts.node(p_keep).time
    for p in p_moves:
        t_move = merged_ts.node(p).time
        delta = t_move - t_keep
        if not (preserve_branch_lengths_when_possible and delta >= 0 and roots_to_move[p]):
            continue
        nodes_to_shift = set()
        for r in roots_to_move[p]:
            q = deque([r])
            while q:
                u = q.popleft()
                if u in nodes_to_shift:
                    continue
                nodes_to_shift.add(u)
                for v in children_of.get(u, ()):
                    q.append(v)
        # Shift node times forward by delta (subtract delta)
        for u in nodes_to_shift:
            tables.nodes.time[u] = tables.nodes.time[u] - delta
        # Optionally shift mutation times
        if shift_mutation_times and tables.mutations.num_rows > 0:
            for mid in range(tables.mutations.num_rows):
                n = tables.mutations.node[mid]
                if n in nodes_to_shift:
                    mt = tables.mutations.time[mid]
                    if mt is not None and not (isinstance(mt, float) and math.isnan(mt)):
                        tables.mutations.time[mid] = mt - delta
        # Ensure moved children are younger than new parent
        for r in roots_to_move[p]:
            if not (tables.nodes.time[r] < t_keep - epsilon):
                raise ValueError(
                    f"Shifted child {r} (time={tables.nodes.time[r]}) not younger than new parent {p_keep} (time={t_keep})"
                )
    
    # Rebuild edge table: reparent and drop specific branches
    to_add_edges = []
    kept_edges = []
    for e in tables.edges:
        if e.parent in p_moves:
            mv = move_by_parent[e.parent]
            if e.child == mv:
                # Drop p_move -> move_leaf
                continue
            if e.child in roots_to_move[e.parent]:
                # Reattach under p_keep
                to_add_edges.append((e.left, e.right, p_keep, e.child, e.metadata))
                continue
            kept_edges.append(e)
            continue
        # Drop root -> p_move edges
        if e.parent == root and e.child in p_moves:
            continue
        kept_edges.append(e)
    
    tables.edges.clear()
    for e in kept_edges:
        tables.edges.add_row(left=e.left, right=e.right, parent=e.parent, child=e.child, metadata=e.metadata)
    for left, right, parent, child, md in to_add_edges:
        if not (tables.nodes.time[parent] > tables.nodes.time[child] + epsilon):
            raise ValueError(
                f"Invalid parent-child times after reparenting: parent {parent} (time={tables.nodes.time[parent]}) "
                f"<= child {child} (time={tables.nodes.time[child]})"
            )
        tables.edges.add_row(left=left, right=right, parent=parent, child=child, metadata=md)
    
    # Update keep_leaf metadata: add abundances from all move_leafs
    keep_md = dict(merged_ts.node(keep_leaf).metadata or {})
    total_add = 0
    for u in move_leafs:
        md = merged_ts.node(u).metadata or {}
        try:
            total_add += int(md.get('abundance', 0))
        except Exception:
            pass
    try:
        base = int(keep_md.get('abundance', 0))
    except Exception:
        base = 0
    keep_md['abundance'] = base + total_add  # stays integer per schema
    
    # Optionally clear sample flags for move_leafs
    flags = tables.nodes.flags.copy()
    if remove_move_leaf_sample_flag:
        for u in move_leafs:
            flags[u] = flags[u] & int(not tskit.NODE_IS_SAMPLE)
        # If you labeled only leaves, you may want to drop labels for move_leafs
        if label_mode == "leaves":
            for u in move_leafs:
                if u in node_labels:
                    del node_labels[u]
    
    # Rebuild the NodeTable row-by-row using metadata dicts (no manual buffers)
    node_tbl = tables.nodes
    schema = node_tbl.metadata_schema
    times = node_tbl.time.copy()       # includes any shifts
    pops = node_tbl.population.copy()
    indivs = node_tbl.individual.copy()
    original_num_rows = node_tbl.num_rows
    
    node_tbl.clear()
    node_tbl.metadata_schema = schema
    for i in range(original_num_rows):
        if i == keep_leaf:
            md_dict = keep_md
        else:
            md_dict = merged_ts.node(i).metadata or {}
        node_tbl.add_row(
            flags=flags[i],
            time=times[i],
            population=pops[i],
            individual=indivs[i],
            metadata=md_dict,
        )
    
    tables.sort()
    final_ts = tables.tree_sequence()
    return(final_ts)


def merge_by_root_and_combine_ancestrals(
ts_list,
tree_index=0,
prefixes=None,
sequence_length=1.0,
root_time=None,
label_mode="all",
root_metadata=None,
site_ancestral_state=None,
site_metadata=None,
preserve_branch_lengths_when_possible=True,
shift_mutation_times=True,
remove_move_leaf_sample_flag=True,
epsilon=1e-9,
):
    if len(ts_list) == 1:
        return ts_list[0]
    
    if prefixes is None:
        prefixes = [f"{i}_" for i in range(len(ts_list))]
    if len(prefixes) != len(ts_list):
        raise ValueError("prefixes length must match number of inputs")
    
    def get_tree(ts, idx):
        if idx < 0 or idx >= ts.num_trees:
            raise IndexError(f"tree_index {idx} out of range for input with {ts.num_trees} trees")
        lefts = list(ts.breakpoints())
        return ts.at(lefts[idx])
    
    # Collect trees and root times
    trees, roots, root_times = [], [], []
    for i, ts in enumerate(ts_list):
        t = get_tree(ts, tree_index)
        rts = t.roots if hasattr(t, "roots") else ([t.root] if t.root != tskit.NULL else [])
        if len(rts) != 1:
            raise ValueError(f"Input {i} has {len(rts)} roots; expected 1")
        r = rts[0]
        trees.append(t)
        roots.append(r)
        root_times.append(ts.node(r).time)
    
    # Common root time and time shifts
    T = max(root_times) if root_time is None else float(root_time)
    shifts = [T - rt for rt in root_times]
    
    tables = tskit.TableCollection(sequence_length=float(sequence_length))
    
    # Metadata schemas
    tables.edges.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {"type": "object", "properties": {"distance": {"type": "integer"}}, "additionalProperties": False},
    })
    tables.nodes.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {
            "type": "object",
            "properties": {"fitness": {"type": "float"}, "abundance": {"type": "integer"}, "mutation_count": {"type": "integer"}},
            "additionalProperties": False,
        },
    })
    tables.mutations.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {"type": "object", "properties": {"s": {"type": "float"}}, "additionalProperties": False},
    })
    tables.sites.metadata_schema = tskit.MetadataSchema({
        "codec": "json",
        "schema": {"type": "object", "properties": {}, "additionalProperties": False},
    })
    
    # Single site (position 0)
    if site_ancestral_state is None:
        site_ancestral_state = ""
        for ts in ts_list:
            if ts.num_sites > 0:
                site_ancestral_state = ts.tables.sites[0].ancestral_state
                break
    if site_metadata is None:
        for ts in ts_list:
            if ts.num_sites > 0:
                md = ts.tables.sites[0].metadata
                if md not in (None, {}, b"", "", []):
                    site_metadata = md
                    break
    if site_metadata is None:
        site_metadata = {}
    site0_id = tables.sites.add_row(position=0.0, ancestral_state=site_ancestral_state, metadata=site_metadata)
    assert site0_id == 0
    
    # Synthetic root node
    root_out_md = {"fitness": np.nan,
                "abundance": 0,
                "mutation_count": 0} if root_metadata is None else dict(root_metadata)
    root_out = tables.nodes.add_row(time=T, flags=0, population=tskit.NULL, individual=tskit.NULL, metadata=root_out_md)
    
    # Map input node IDs to output node IDs; build labels
    node_maps = []
    node_labels = {}
    for i, ts in enumerate(ts_list):
        m = np.full(ts.num_nodes, -1, dtype=np.int32)
        prefix = prefixes[i]
        shift = shifts[i]
        is_sample = np.zeros(ts.num_nodes, dtype=bool)
        is_sample[ts.samples()] = True
        for u in range(ts.num_nodes):
            n = ts.node(u)
            new_id = tables.nodes.add_row(
                time=n.time + shift,
                flags=n.flags,
                population=tskit.NULL,
                individual=tskit.NULL,
                metadata=n.metadata,
            )
            m[u] = new_id
            if label_mode == "all" or (label_mode == "leaves" and is_sample[u]):
                node_labels[new_id] = f"{prefix}{u}"
        node_maps.append(m)
    
    # Add edges over full interval; reparent input roots to synthetic root
    SL = float(sequence_length)
    for i, (ts, tree, r_in) in enumerate(zip(ts_list, trees, roots)):
        m = node_maps[i]
        for u in tree.nodes():
            p = tree.parent(u)
            if p == tskit.NULL:
                continue
            e_id = getattr(tree, "edge", lambda _: tskit.NULL)(u)
            md = ts.tables.edges[e_id].metadata if e_id != tskit.NULL else None
            parent_out = root_out if p == r_in else m[p]
            child_out = m[u]
            tables.edges.add_row(left=0.0, right=SL, parent=parent_out, child=child_out, metadata=md)
    
    # Copy mutations
    for i, ts in enumerate(ts_list):
        m_map = node_maps[i]
        has_mut_time = hasattr(ts.tables.mutations, "time") and len(ts.tables.mutations.time) == ts.num_mutations
        shift = shifts[i]
        mut_id_map = np.full(ts.num_mutations, -1, dtype=np.int32)
        pending = list(range(ts.num_mutations))
        while pending:
            progressed = False
            next_pending = []
            for old_mid in pending:
                m = ts.tables.mutations[old_mid]
                parent_ok = (m.parent == tskit.NULL) or (mut_id_map[m.parent] != -1)
                if parent_ok:
                    new_parent = tskit.NULL if m.parent == tskit.NULL else mut_id_map[m.parent]
                    new_mid = tables.mutations.add_row(
                        site=0,
                        node=m_map[m.node],
                        derived_state=m.derived_state,
                        parent=new_parent,
                        time=(m.time + shift) if has_mut_time and m.time is not None else None,
                        metadata=m.metadata,
                    )
                    mut_id_map[old_mid] = new_mid
                    progressed = True
                else:
                    next_pending.append(old_mid)
            if not progressed:
                raise ValueError("Could not resolve mutation.parent references; cycle or missing parent?")
            pending = next_pending
    
    # Build temporary merged TreeSequence
    tables.sort()
    merged_ts = tables.tree_sequence()
    
    # Identify keep_leaf and move_leafs
    tree = merged_ts.first()
    root = tree.root
    if root == tskit.NULL:
        raise ValueError("Merged tree has no root")
    
    def is_terminal_sample(u):
        return ((merged_ts.node(u).flags & tskit.NODE_IS_SAMPLE) != 0) and (tree.num_children(u) == 0)
    
    def md_ok(u):
        md = merged_ts.node(u).metadata or {}
        return (md.get('abundance', 0) > 0) and (md.get('mutation_count', 0) == 0)
    
    leaves_in_order = list(tree.leaves(root))
    candidates = [u for u in leaves_in_order if is_terminal_sample(u) and md_ok(u)]
    if len(candidates) <= 1:
        return merged_ts
    
    keep_leaf = candidates[0]
    move_leafs = candidates[1:]
    
    # Allow keep_leaf to be either direct child of root or under an internal parent
    target_parent = tree.parent(keep_leaf)
    if target_parent == tskit.NULL:
        # Degenerate case: keep_leaf is itself a root; nothing sensible to merge
        return merged_ts
    
    # Partition move_leafs into two groups
    internal_p_moves = set()  # parents that are children of the root
    move_by_parent = {}
    root_children_moves = set()  # move_leafs whose parent is the root
    # print(merged_ts.draw_text())
    for u in move_leafs:
        p = tree.parent(u)
        if p == tskit.NULL:
            continue
        if p == root:
            root_children_moves.add(u)
        elif tree.parent(p) == root:
            internal_p_moves.add(p)
            move_by_parent[p] = u
        else:
            # Not relaxing beyond one level; keep behavior unchanged for deeper parents
            raise ValueError(f"move_leaf {u} must be under the root or under a node whose parent is the root")
    
    # Dump tables for editing
    tables = merged_ts.dump_tables()
    
    # children map
    children_of = defaultdict(set)
    for e in tables.edges:
        children_of[e.parent].add(e.child)
    
    # For internal parents, determine children to reattach under target_parent
    roots_to_move = {}
    for p in internal_p_moves:
        mv = move_by_parent[p]
        roots_to_move[p] = [c for c in children_of[p] if c != mv]
    
    # Time shifting for internal parents relative to target_parent
    t_target = merged_ts.node(target_parent).time
    for p in internal_p_moves:
        t_move = merged_ts.node(p).time
        delta = t_move - t_target
        if not (preserve_branch_lengths_when_possible and delta >= 0 and roots_to_move[p]):
            continue
        nodes_to_shift = set()
        for r in roots_to_move[p]:
            q = deque([r])
            while q:
                u = q.popleft()
                if u in nodes_to_shift:
                    continue
                nodes_to_shift.add(u)
                for v in children_of.get(u, ()):
                    q.append(v)
        # Shift node times forward by delta
        for u in nodes_to_shift:
            tables.nodes.time[u] = tables.nodes.time[u] - delta
        # Optionally shift mutation times
        if shift_mutation_times and tables.mutations.num_rows > 0:
            for mid in range(tables.mutations.num_rows):
                n = tables.mutations.node[mid]
                if n in nodes_to_shift:
                    mt = tables.mutations.time[mid]
                    if mt is not None and not (isinstance(mt, float) and math.isnan(mt)):
                        tables.mutations.time[mid] = mt - delta
        # Ensure moved children younger than target_parent
        for r in roots_to_move[p]:
            if not (tables.nodes.time[r] < t_target - epsilon):
                raise ValueError(
                    f"Shifted child {r} (time={tables.nodes.time[r]}) not younger than new parent {target_parent} (time={t_target})"
                )
    
    # Rebuild edges
    to_add_edges = []
    kept_edges = []
    for e in tables.edges:
        # Internal parent case: drop p_move->move_leaf, reattach siblings, drop root->p_move
        if e.parent in internal_p_moves:
            mv = move_by_parent[e.parent]
            if e.child == mv:
                continue  # drop p_move -> move_leaf
            if e.child in roots_to_move[e.parent]:
                # Reattach under target_parent (could be the root)
                to_add_edges.append((e.left, e.right, target_parent, e.child, e.metadata))
                continue
            kept_edges.append(e)
            continue
        # Drop root -> p_move edges for internal parents
        if e.parent == root and e.child in internal_p_moves:
            continue
        # Root-children move case: drop root -> move_leaf
        if e.parent == root and e.child in root_children_moves:
            continue
        kept_edges.append(e)
    
    tables.edges.clear()
    for e in kept_edges:
        tables.edges.add_row(left=e.left, right=e.right, parent=e.parent, child=e.child, metadata=e.metadata)
    for left, right, parent, child, md in to_add_edges:
        if not (tables.nodes.time[parent] > tables.nodes.time[child] + epsilon):
            raise ValueError(
                f"Invalid parent-child times after reparenting: parent {parent} (time={tables.nodes.time[parent]}) "
                f"<= child {child} (time={tables.nodes.time[child]})"
            )
        tables.edges.add_row(left=left, right=right, parent=parent, child=child, metadata=md)
    
    # Update keep_leaf abundance
    keep_md = dict(merged_ts.node(keep_leaf).metadata or {})
    total_add = 0
    for u in move_leafs:
        md = merged_ts.node(u).metadata or {}
        try:
            total_add += int(md.get('abundance', 0))
        except Exception:
            pass
    try:
        base = int(keep_md.get('abundance', 0))
    except Exception:
        base = 0
    keep_md['abundance'] = base + total_add
    
    # Optionally clear sample flags for all move_leafs
    flags = tables.nodes.flags.copy()
    # if remove_move_leaf_sample_flag:
    #     for u in move_leafs:
    #         flags[u] &= ~tskit.NODE_IS_SAMPLE
    
    idx = np.array(move_leafs, dtype=np.int64)
    flags[idx] &= ((~tskit.NODE_IS_SAMPLE) & 0xFFFFFFFF)
    
    # Rebuild NodeTable row-by-row using metadata dicts
    node_tbl = tables.nodes
    schema = node_tbl.metadata_schema
    times = node_tbl.time.copy()
    pops = node_tbl.population.copy()
    indivs = node_tbl.individual.copy()
    original_num_rows = node_tbl.num_rows
    
    node_tbl.clear()
    node_tbl.metadata_schema = schema
    for i in range(original_num_rows):
        if i == keep_leaf:
            md_dict = keep_md
        else:
            md_dict = merged_ts.node(i).metadata or {}
        node_tbl.add_row(
            flags=flags[i],
            time=times[i],
            population=pops[i],
            individual=indivs[i],
            metadata=md_dict,
        )
    
    tables.sort()
    final_ts = tables.tree_sequence()
    return final_ts