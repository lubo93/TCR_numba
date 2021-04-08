# info: try to use the numba_test_2.py created sparse file to make a network

import numpy as np
#import networkx as nx
import matplotlib.pyplot as plt
import time

from graph_tool.all import *
import graphics

from matplotlib import rcParams

import pandas as pd 

import seaborn as sns

def convert_edges(name="data/sparse.txt", len_x=10, single_sidelength=4*10**3):
    """
    Info: open file with name "name" and convert the contained single number indices
        to two-number indices
    Args: name: name of the file
        len_x: number of blocks in a single row
        single_sidelength: sidelength of a single block
    Returns: G: assembled graph
    """
    xs = list()
    ys = list()

    sidelength = single_sidelength*len_x
    
    data = None
    
    with open(name, "r") as f:
        lines = f.readlines()
        #data = np.zeros([len(lines), 2])
        data = list()
        for i in range(len(lines)):
            #data = np.array([line.strip().split() for line in f],float)
            var = int(lines[i])
            #print("\n var: ", var, "; var.strip(): ", var.strip(), "; var.strip().split(): ", var.strip().split())
            #var = int(var.strip())
            val1 = var%sidelength
            val2 = int(var/sidelength)
            #data.append([val1, val2])
            if val1 != val2:
                #data[i][0] = val1
                #data[i][1] = val2
                data.append([val1, val2])
                #print("\n unequal val1: ", val1, "; val2: ", val2)
            else:
                print("\n data i: ", i, "; equal val1: ", val1, "; val2: ", val2)
        print("\n L: ", len(lines))
        #lines = np.loadtxt('sparse.txt')
    return data

def assemble_graph(data, cut_value = -1):
    """
    Info: assemble graph from two-number edge list "data"
    Args: data: two-number edgelist
        cut_value: int number: In case you only want to use the first 
            N edges of the file, you can set input some cut_value=N, otherwise 
            it should usually be set to -1.
    Returns: G: assembled graph
    """

    # info: make a graph and 
    G = Graph()
    ug = Graph(directed = False)
    # info: adding nodes, in order to first create all the nodes for the network
    t0 = time.time()
    #G.add_nodes_from(list(range(len(data))))
    #for i )):
    t1 = time.time()
    
    # info: adding the edges
    #G.add_edges_from(data)
    vals = list()
    dic = {}
    len_data = len(data)
    vals_edge_list = list()

    # TODO: D the analysis for higher index
    if cut_value == -1:
        len_val = len(data)
    else:
        len_val = min(len(data), cut_value)
    
    for i in range(len_val):#(len(data)):#range(len(data)):
        data_i0 = data[i][0]
        data_i1 = data[i][1]
        
        try: 
            idx_1 = dic[data_i0]
        except: 
            idx_1 = len(list(dic.keys()))
            dic[data_i0] = idx_1

        try:
            idx_2 = dic[data_i1]
            vals_edge_list.append([idx_2])
        except:
            idx_2 = len(list(dic.keys()))
            dic[data_i1] = idx_2

        vals_edge_list.append([idx_1, idx_2])
        """
        if data_i0 in vals: 
            pass
        else: 
            vals.append(data_i0)
        
        if data_i1 in vals:
            pass
        else:
            vals.append(data_i1)
        vals_edge_list.append([vals.index(data_i0), vals.index(data_i1)])
        """
        #G.add_edge(vals.index(data_i0), vals.index(data_i1))
        
        #G.add_edge(data[i][0], data[i][1])
        if i%1000 == 0:
            print("\n i: ", i, " / ", len_data, "; data[i][0]: ", data[i][0], "; data[i][1]: ", data[i][1], "; len(list(G.vertices())): ",  len(list(G.vertices())), "; dt: ", (time.time() - t1))
        if i%10**4 == 0:
            print("\n dt: ", (time.time() - t1))
    print("\n making graph ... ")
    print("\n len(vals_edge_list): ", len(vals_edge_list))
    G.add_edge_list(vals_edge_list)
    print("\n graph finished")
    t2 = time.time()
        
    for v in G.vertices():
        neighbors = v.all_neighbors()
        #print("\n v: ", [G.vertex_index[neighbor] for neighbor in neighbors])
    #print("\n len(data): ", len(data)) 
    #print("\n len(list(G.vertices())): ",  len(list(G.vertices())))
    """   
    # info: if the graph is empty, then add at least one node,
    #     which is called the "ZeroNode"
    if len(G.nodes()) == 0:
        G = nx.Graph()
        G.add_node("ZeroNode")
    """
    return G

def degree_distribution(G):
    """
    Info: calculate and return the degree distribution for 
        the graph G
    Args: G: graph
    Returns: deg_freqs: dictionary with degree values 
        as keys and their frequencies as values
    """
    
    deg_freqs = {}

    for v in G.vertices():
        deg = G.get_total_degrees([v])
        deg = deg[0]
        try:
            deg_freqs[deg] += 1
        except: 
            deg_freqs[deg] = 1
    return deg_freqs

def clique_distribution(g):
    """
    Info: calculate and return clique distribution for 
        the graph g
    Args: g: graph
    Returns: xs: array of cliquesized
        ys: corresponding array of frequencies for each 
            clique size
    """

    cliques = list(max_cliques(g))
    sz = [len(list(clique)) for clique in cliques]

    dic = {}

    for i in range(len(sz)):
        l = sz[i]
        try: 
            dic[l] += 1
        except: 
            dic[l] = 1
    xs = dic.keys()
    ys = dic.values()
    
    return xs, ys

def neighbor_degree(g):
    """
    Info: calculate and return the neighbor degree distribution
        for graph g
    Args: g: Graph
    Returns: xs: array of degree values
        ys: corresponding average neighbor degrees
    """
    # perh. TODO: check that the function really does, what we say in 
    #     the description

    dic = {}
    
    len_total = len(list(g.vertices()))

    print("\n long ?")
    for v in g.vertices():
        neighbors = list(v.all_neighbors())
        k_self = len(neighbors)
        neighbor_idxs = [g.vertex_index[neighbor] for neighbor in neighbors]
        k_neighbors = [g.get_total_degrees([w])[0] for w in neighbor_idxs]
        k_mean = np.mean(k_neighbors)
        try:
            dic[k_self].append(k_mean)
        except:
            dic[k_self] = [k_mean]
    print("\n long done")

    print("\n len(dic.keys()): ", len(dic.keys()))
    for key in dic.keys():
        dic[key] = np.mean(dic[key])

    print("\n keys done")
    xs = dic.keys()
    ys = dic.values()

    return xs, ys

def robustness(g):
    """
    Info: perform the robustness anaplysis for a 
        graph  g and return the results
    Args: g: Graph
    Returns: xs: array remaining node counts
        ys: corresponding maximum  cluster sizes
    """
    
    vertices = sorted([v for v in g.vertices()])
    sizes, comp = vertex_percolation(g, vertices)
    
    xs = list(range(len(sizes)))
    ys = sizes

    return xs, ys
 
def plot(xs, ys, log=True, name=None):
    """
    Info: plot the array ys over the array xs with 
        some standard settings and save it
    Args: xs: array of x-values
        ys: array of y-values
        log: bool: True: log-log-plot; False: no log-plot
        name: name of the file to save to
    Returns: -
    """
    
    fig_size, params = graphics.design(factor=1)
    rcParams.update({'figure.figsize': fig_size})

    fig_1, ax_1 = plt.subplots()

    ax_1.plot(xs, ys, '.')
    if log == True:
        ax_1.set_xscale('log')
        ax_1.set_yscale('log')
    if name!=None:
        fig_1.tight_layout()
        fig_1.savefig(name, dpi=300, bbox_inches="tight")
    plt.show()

def pathogen_analysis(g, idx_max):
    """
    Info: perform the network analysis for the graph g
    Args: g: the graph (supposed to originate from virus 
            included edge list)
        idx_max: number of virus cognate CDR3 sequences
    Returns: res: an object, containing all the results
    """

    res = {}
    
    vp, ep = betweenness(g)
    vp_list = list(vp)
    res["betweenness"] = vp_list

    dic = {}
    for i in range(idx_max):
        try: 
            dic[vp_list[i]] += 1
        except: 
            dic[vp_list[i]] = 0
    res["betweenness_dist"] = dic

    return res

def main(): 
    """
    Info: Do all the graph analysis
    Args: -
    Returns: -
    """

    """
    g = assemble_graph(, cut_value=4*10**4)
    deg_freqs = degree_distribution(g)
    xs = list()
    ys = list()
    
    # info: make degree distribution
    for key in deg_freqs.keys():
        xs.append(key)
        ys.append(deg_freqs[key])
    plot(xs, ys)
    
    print("\n cliques ...")
    # info: make clique size distribution
    
    xs, ys = clique_distribution(g)
    name = "new_plots/cliques.png"
    plot(xs, ys, name=name)

    print("\n neighbor_degree ...")
    # info: neighbor degree: 
    xs, ys = neighbor_degree(g)
    name = "new_plots/neighbors.png"
    plot(xs, ys, name=name)

    print("\n robustness ...")
    xs, ys = robustness(g)
    name = "new_plots/robustness.png"
    plot(xs, ys, log=False, name=name)
    """
    print("\n pathogens ...")
    # info: pathogens
    pathogens = ["Influenza", "HIV", "Alzheimer", "Parkinson", "Tuberculosis", "COVID-19"]
    cols = {"Influenza": (1.0, 0.0, 0.0, 1.0), 
            "HIV": (0.0, 1.0, 0.0, 1.0), 
            "Alzheimer": (0.0, 0.0, 1.0, 1.0),
            "Parkinson": (0.0, 0.0, 0.0, 1.0), 
            "Tuberculosis": (0.5, 0.5, 0.5, 1.0), 
            "COVID-19": (0.5, 0.5, 0.0, 1.0)}
    ress = {}
    
    ress_pathogen = list()
    ress_val = list()
    
    fig, ax = plt.subplots(nrows=2, ncols=2)
    ax_00 = ax[0][0]
    ax_01 = ax[0][1]
    ax_10 = ax[1][0]
    ax_11 = ax[1][1]
    
    for pathogen in pathogens:
        name = "data/" + pathogen + "_small.txt"
        data = convert_edges(name=name, len_x=10, single_sidelength=4*10**3)
        g = assemble_graph(data)
        idx_max = 100
        res_raw = pathogen_analysis(g, idx_max)
        #TODO: load idx_max value
        res = res_raw["betweenness"][:idx_max]
        print("\n pathogen: ", pathogen, "; res: ", res)
        for i in range(len(res)):
            ress_pathogen.append(pathogen)
            ress_val.append(res[i])
        dic = res_raw["betweenness_dist"]
        keys = list(dic.keys())
        vals = list(dic.values())
        
        ax_10.plot(keys, vals, '.', color=cols[pathogen])
    #res[pathogen] = res["betweenness"] 
    plt.show()   
    #xs, ys = viruses(g)
    ress["pathogen"] = ress_pathogen
    ress["val"] = ress_val
    data = pd.DataFrame(ress)
    sns.stripplot(data=data, x="pathogen", y="val", ax=ax_00)

    #data = pd.DataFrame(ress)
    #ax_01.plot(ress)
    """
    g = assemble_graph("sparse_virus.txt")
    """
    # TODO: some virus graph analysis

if __name__ == "__main__":
    print("\n name: ", name)
    main()
