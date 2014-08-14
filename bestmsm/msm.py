""" 
================
 MSM 
================
"""

import sys
import types
import numpy as np
import networkx as nx


class MSM:
    """
    A class for constructing the MSM

    Parameters
    ----------
    trajectories : list of str
        Set of trajectories used for the construction.

    lag : float
        Lag time for building the MSM.

    Attributes
    ----------
    keys : list of str
        Names of states.

    count : np array
        Transition count matrix.

    keep_states : list of str
        Names of states after removing not strongly connected sets.
        
    """
    def __init__(self, trajectories):
        self.keep_states = None
        self.keep_keys = None
        # merge state keys from all trajectories
        self.keys = self.merge_trajs(trajectories)
        self.data = trajectories
#        # find largest strongly connected sets
#        self.keep_states, self.keep_keys = self.check_connect()
        # calculate transition matrix
#        self.T = self.calc_trans()
#        self.K = self.calc_rate()
        # calculate eigenvalues and eigenvectors
#        self.calc_eigs()

    def merge_trajs(self, trajectories):
        """ Merge all trajectories into a consistent set. """
        new_keys = []
        for traj in trajectories:
            new_keys += filter(lambda x: x not in new_keys, traj.keys)
        return new_keys
    
    def calc_count(self, lagt=None):
        """ Calculate transition count matrix. """
        # initialize count matrix 
        keys = self.keys
        nkeys = len(keys)
        count = np.zeros([nkeys,nkeys], int)

        # loop over trajectories
        for traj in self.data:
            raw = traj.states
            dt = traj.dt
            lag = int(lagt/dt) # lag time in frames
            nraw = len(raw)
            for i in range(0, nraw-lag, lag):
                j = i + lag
                state_i = raw[i]
                state_j = raw[j]
                if state_i in keys:
                    idx_i = keys.index(state_i)
                if state_j in keys:
                    idx_j = keys.index(state_j)
                try:
                    count[idx_j][idx_i] += 1
                except IndexError:
                    pass
        return count
    
    def check_connect(self):
        """ Check connectivity of rate matrix. """
        D = nx.DiGraph(self.count)
        keep_states = sorted(nx.strongly_connected_components(D)[0])
        keep_keys = map(lambda x: self.keys[x], self.keep_states)
        return keep_states, keep_keys

    def calc_trans(self, count=None, lagt=None):
        """ Calculate transition matrix """
        if not isinstance(count, list):
            count = self.calc_count(lagt)
            print count
        keep_states = self.keep_states
        nkeep = len(keep_states)
        T = np.zeros([nkeep,nkeep], float)
        for i in range(nkeep):
            ni = reduce(lambda x, y: x + y,
                        map(lambda x: count[keep_states[x]][keep_states[i]],
                            range(nkeep)))
            for j in range(nkeep):
                T[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)
        return T
    
    def calc_rate(self, lagt):
        """ Calculate rate matrix using a Taylor series as described in De Sancho et al
        J Chem Theor Comput (2013)"""
        nkeep = len(keep_states)
        K = self.T/lagt
        for i in range(nkeep):
            K[i][i] = -(np.sum(K[:i,i]) + np.sum(K[i+1:,i]))
        return K

    def calc_eigs(self):
        """ Calculate eigenvalues and eigenvectors"""
        #evalsK,rvecsK = np.linalg.eig(K)
        self.evalsK, self.lvecsK, self.rvecsK = scipyla.eig(self.K, left=True)
        self.evalsT, self.lvecsT, self.rvecsT = scipyla.eig(self.T, left=True)
        # from K
        elistK = []
        for i in range(nkeep):
            elistK.append([i,np.real(evalsK[i])])
        elistK.sort(msm_lib.esort)
        ieqK,eK = elistK[0]
        peqK_sum = reduce(lambda x,y: x + y, map(lambda x: rvecsK[x,ieqK], 
            range(nkeep)))
        peqK = rvecsK[:,ieqK]/peqK_sum
        # from T
        elistT = []
        for i in range(nkeep):
            elistT.append([i,np.real(evalsT[i])])
        elistT.sort(msm_lib.esort)
        ieqT,eT = elistT[0]
        peqT_sum = reduce(lambda x,y: x + y, map(lambda x: rvecsT[x,ieqT],
             range(nkeep)))
        peqT = rvecsT[:,ieqT]/peqT_sum
        return peqT, rvecsT
