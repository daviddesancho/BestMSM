import sys
import math
import itertools
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# Much indebted to Randal Olson (http://www.randalolson.com/)
# These are the "Tableau 20" colors as RGB.  
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)

figsize = (12,14)


def getridofticks():
    """ Ticking only at bottom and left of plot"""
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def write_dot(J, nodeweight=None, rank=None, out="out.dot"):
    """ Function for printing a graph in dot format
    
    Parameters
    ----------
    J : array
        A directed graph in the form of a matrix, with values corresponding to weights
    nodeweight : array
        Weights for the nodes.
    rank : list, array
        The rank order for display of nodes.
    out : string
        Filename for output.
    
    """
    print "\n Writing graph in dot format..."
    fout = open(out, "w")
    fout.write("strict digraph G {\n")

    D = nx.DiGraph(J.transpose())
    nd = len(D.nodes())

    # define weights of nodes
    try:        
        lw = np.log(nodeweight)
        #wmin = np.min(lw)
        #wmax = np.max(lw)
        #norm = wmax - wmin
        #weight = [(x - wmin)/norm + 0.5 for x in lw]
        weight = 0.2*(lw - np.min(lw) -np.max(lw))
    except AttributeError:
        weight = np.ones(nd)

    elems = zip(*D.edges())[0] + zip(*D.edges())[1]
    
    for n in D.nodes():
        if n in elems:
            fout.write("%i [shape=circle,width=%f];\n"%(n, weight[n]))

    # define rank of nodes
    try:
        for u,v in itertools.product(D.nodes(),D.nodes()):
            if u  < v:
                if u in elems and v in elems:
                    if rank[u] == rank[v]:
                        fout.write ("{rank=same; %i; %i;}\n"%(u,v))
    except TypeError:
        pass

    # define weights of edges
    weight = np.zeros((nd, nd), float)
    try:
        wmin = np.min(J[J>0])
        for u, v in D.edges():
            D[u][v]['weight'] = J[v,u]/wmin
    except TypeError:
        for u, v in D.edges():
            D[u][v]['weight'] = 1. 

    for u,v,d in D.edges_iter(data=True):
        if u != v and not math.isnan(d['weight']):
            fout.write("%i -> %i  [penwidth=%f,color=black];\n"%(u, v, np.log(d['weight'])+1))

    fout.write("}")

def plot_evals(vals, log=False):
    """ Plot eigenvalues """
    fig = plt.figure()
    ax = plt.add_subplots(1,1,1)
    getridofticks()
    if not log:
        plt.plot(vals)
    else:
        plt.semilogy(vals)
    ax.set_xlabel(r'$Eigenvalue$')
    ax.set_ylabel(r'$\tau$')
    plt.show()
