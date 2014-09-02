import numpy as np
import matplotlib.pyplot as plt

def write_dot(D, nodeweight=None, edgeweight=None, out="out.dot"):
    """ Function for printing a graph in dot format
    
    Parameters
    ----------
    D : array
        A directed graph.
    nodeweight : array
        Weights for the nodes.
    edgeweight : array 
        Weights for the edges.
    out : string
        Filename for output.
    
    """
    print "\n Writing graph in dot format..."
    fout = open(out, "w")
    fout.write("strict digraph G {\n")
    print D
    nd = len(D.nodes())
    # define weights of nodes
    try:        
        lw = np.log(nodeweight)
        wmin = np.min(lw)
        wmax = np.max(lw)
        norm = wmax - wmin
        weight = [(x - wmin)/norm + 0.5 for x in lw]
    except AttributeError:
        weight = np.ones(nd)
    for n in D.nodes():
        fout.write("%i [shape=circle,width=%f];\n"%(n, weight[n]))

    # define weights of edges
    weight = np.zeros((nd, nd), float)
    try:
        lw = [np.log(x) for row in edgeweight for x in row if x > 0]
        wmin = np.min(lw)
        wmax = np.max(lw)
        norm = wmax - wmin
        for u, v in D.edges():
            weight[u,v] = 2*((np.log(edgeweight[u,v]) - wmin)/norm + 0.5)
    except TypeError:
        weight = np.ones((nd, nd), float)
#    #   if pfold[u] == pfold[v]:
#    #       out.write ("{rank=same; %i; %i;}\n"%(u,v))
#

    for (u,v) in D.edges(): 
       fout.write("%i -> %i  [penwidth=%f,color=black];\n"%(u, v, weight[u,v]))

    fout.write("}")
