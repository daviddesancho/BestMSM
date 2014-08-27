import sys
from msm import MSM
import trajectory as traj
import pcca_lib
import copy

""" 
================
 PCCA 
================
"""

class PCCA(MSM):
    """
    A class for doing PCCA clustering

    Parameters
    ----------
    N : int
        The desired number of clusters.

    Attributes
    ----------
    keys : dict
        A dictionary containing the clusters formed.

    """

    def __init__(self, parent, N=2, robust=True, optim=False):
        keys = range(N)
        self.parent = parent
        self.macros = self.eigen_group(N=N)
        self.N = N
        self.mappedtraj = self.map_trajectory()
        lagt = parent.lagt
        super(PCCA, self).__init__(self.mappedtraj, keys=keys, lagt=lagt)

    def eigen_group(self, N=2, method="robust"):
        """ Splits microstates into macrostates 
        
        Parameters
        ----------
        N : int
            Number of clusters.

        method : str
            The method used for clustering.

        Returns
        -------
        macros : dict
            A dictionary with the membership to macrostates.

        """

        # generate eigenvectors in case the MSM does not have them
        if not hasattr(self.parent, 'lvecsT'):
           lagt = self.parent.lagt
           tauT, peqT, self.parent.rvecsT, self.parent.lvecsT = \
                   self.parent.calc_eigs(lagt=lagt, evecs=True)
        lvecs = self.parent.lvecsT

        # split in desired number of macrostates
        macros = {}
        macros[0] = self.parent.keep_states 
        for n in range(1, N):
            if method is "robust":
                macro_new, vals = pcca_lib.split_sigma(macros, lvecs[:,n])
            elif method is "sign":
                macro_new, vals = pcca_lib.split_sign(macros, lvecs[:,n])
            macros = copy.deepcopy(macro_new)
        print "\n Initial membership of microstates to macrostates:"
        if len(self.parent.keep_keys) < 100:
            for k,v in macros.iteritems():
                print k, [self.parent.keep_keys[x] for x in v]
        else:
            for k,v in macros.iteritems():
                print k,":", len(v)
        return macros

    def map_trajectory(self):
        mappedtraj = []
        for data in self.parent.data:
            mt = traj.TimeSeries()
            mt.dt = data.dt
            mt.time = data.time
            mt.states = []
            for s in data.states:
                mt.states.append([k for k, v in self.macros.iteritems() \
                        if self.parent.keep_keys.index(s) in v][0])
            mappedtraj.append(mt)
        return mappedtraj 
