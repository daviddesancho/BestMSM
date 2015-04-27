import numpy as np
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



