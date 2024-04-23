#!/usr/bin/python
"""Functions for merging multiple graphs."""
# Updated: 12/04/21
import numpy as np
import leidenalg as la


def compute_graphs(mod, plist, keep_all_z=True):
    """Compute graphs along a list of eigenvalue-power values."""
    for power in plist:
        graph_id = f'p{power}'
        if graph_id not in mod.graphs.keys():
            mod.make_graph(graph_id, full_graph_only=True, power=power,
                           keep_all_z=keep_all_z)
    return(['p' + str(x) for x in plist])


def get_intersection_graph_names(mod, graphs):
    """Get the names shared by all graphs."""
    full_nodes = None
    for x in graphs:
        nodes = mod.graphs[x].graph.vs['name']
        if full_nodes is None:
            full_nodes = nodes
        else:
            nodes = set(nodes)
            full_nodes = [x for x in full_nodes if x in nodes]
    return(full_nodes)


def flag_graphlist_nodes(graphlist, min_size=4):
    """Return the list of nodes in a min_size+ sized component in any graph."""
    keep_nodes = []
    for g in graphlist:
        memb = g.clusters().membership  # Connected components
        u, c = np.unique(memb, return_counts=True)  # Sizes of components
        comp = set(u[c < min_size])  # Components to remove
        to_keep = [i for i, x in enumerate(memb) if x not in comp]
        keep_nodes = keep_nodes + to_keep
    return(keep_nodes)


def prune_graph_to_nodes(graph, full_nodes):
    """Remove nodes in a graph not in full_nodes."""
    nodes = graph.vs['name']
    todel = [i for i, x in enumerate(nodes) if x not in full_nodes]
    graph.delete_vertices(todel)
    return(graph)


def clean_graphlist_nodes(graphlist, keep_nodes):
    """Remove nodes in a graph not in keep_nodes."""
    # TODO MERGE with ABOVE
    NV = len(graphlist[0].vs['name'])
    keep_nodes = np.unique(np.array(keep_nodes))
    to_delete = [i for i in np.arange(NV) if i not in keep_nodes]
    for g in graphlist:
        g.delete_vertices(to_delete)  # Remove nodes in components
    return(graphlist)


def make_graphlist(mod, plist, min_size=4, keep_all_z=True):
    """Make and process a list of graphs."""
    # TODO: allow this to work with a variety of values (SVD k, power, z, etc)
    # Compute the graphs on the modules object:
    graphs = compute_graphs(mod, plist, keep_all_z=keep_all_z)
    # Merge gene sets to ensure same sets:
    full_nodes = get_intersection_graph_names(mod, graphs)
    # Prune graphs to same gene set:
    graphlist = []
    for x in graphs:
        g = mod.graphs[x].graph.copy()
        g = prune_graph_to_nodes(g, full_nodes)
        graphlist.append(g)
    # Flag nodes not in large components:
    keep_nodes = flag_graphlist_nodes(graphlist, min_size=min_size)
    # Delete vertices not in keepnodes:
    graphlist = clean_graphlist_nodes(graphlist, keep_nodes)
    return(graphlist, graphs)


def partition_graphlist(graphlist, resolution=3,
                        n_iterations=-1, random_state=1):
    """Partition nodes on multiple graphs using leidenalg multiplex."""
    partition_type = la.RBConfigurationVertexPartition
    partition_kwargs = {}
    partition_kwargs["n_iterations"] = n_iterations
    partition_kwargs["seed"] = random_state
    partition_kwargs["resolution_parameter"] = resolution
    memb, _ = la.find_partition_multiplex(
        graphlist, partition_type, **partition_kwargs)
    return(memb)
