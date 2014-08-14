""" 
================
 MSM 
================
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class MasterMSM:
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
        self.data = trajectories
        self.keys = self.__merge_trajs()
        self.dt = self.__min_dt()
        self.__out()
        self.__chapman_kolmogorov()

    def __merge_trajs(self):
        """ Merge all trajectories into a consistent set"""
        new_keys = []
        for traj in self.data:
            new_keys += filter(lambda x: x not in new_keys, traj.keys)
        return new_keys

    def __min_dt(self):
        """ Find minimum dt in trajectories"""
        return np.min(map(lambda x: x.dt, self.data))

    def __out(self):
        """ Output description to user"""
        print " Building MSM from "
        print map(lambda x: x.filename, self.data)

    def __chapman_kolmogorov(self):
        """ Carry out Chapman-Kolmogorov test"""

        # defining lag times to produce the MSM
        #lagtimes = self.dt*np.array([1, 2, 5, 10, 20, 50, 100, 200, 500])
        lagtimes = self.dt*np.array([1])
        msms = {}
        data = []
        for lagt in lagtimes:
            print "\n Generating MSM at lag time: %g"%lagt
            msms[lagt] = MSM(self.data, self.keys, lagt)
            print "\n    Count matrix:"
            print msms[lagt].count
            print "\n    Transition matrix:"
            print msms[lagt].trans
            data.append(lagt, msms[lagt].tauT[1])
        plt.plot(data[:,0], data[:,1])
        plt.show()


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
    def __init__(self, data, keys=None, lagt=None):
        # merge state keys from all trajectories
        self.keys = keys
        self.data = data
        self.count = self.calc_count(lagt)
        self.keep_states, self.keep_keys = self.check_connect()
        self.trans = self.calc_trans(lagt)
        self.tauT = self.calc_eigs()

    def calc_count(self, lagt=None):
        """ Calculate transition count matrix. """
        # initialize count matrix
        print self.keys
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
        print "\n checking connectivity"
        D = nx.DiGraph(self.count)
        keep_states = sorted(nx.strongly_connected_components(D)[0])
        keep_keys = map(lambda x: self.keys[x], keep_states)
        return keep_states, keep_keys

    def calc_trans(self, lagt=None):
        """ Calculate transition matrix """
        nkeep = len(self.keep_states)
        keep_states = self.keep_states
        count = self.count
        trans = np.zeros([nkeep, nkeep], float)
        for i in range(nkeep):
            ni = reduce(lambda x, y: x + y, map(lambda x: 
                count[keep_states[x]][keep_states[i]], range(nkeep)))
            for j in range(nkeep):
                trans[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)
        return trans
    
    def calc_rate(self, lagt=None):
        """ Calculate rate matrix using a Taylor series 
        as described in De Sancho et al 
        J Chem Theor Comput (2013)
        """
        nkeep = len(self.keep_states)
        rate = self.trans/lagt
        for i in range(nkeep):
            rate[i][i] = -(np.sum(rate[:i,i]) + np.sum(rate[i+1:,i]))
        return rate

    def calc_eigs(self, lagt=None):
        """ Calculate eigenvalues and eigenvectors"""
        #evalsK,rvecsK = np.linalg.eig(K)
        self.evalsK, self.lvecsK, self.rvecsK = \
                   scipyla.eig(self.rate, left=True)
        self.evalsT, self.lvecsT, self.rvecsT = \
                scipyla.eig(self.trans, left=True)
#        # sort modes
        elistK = []
        for i in range(nkeep):
            elistK.append([i,np.real(evalsK[i])])
        elistK.sort(msm_lib.esort)
        elistT = []
        for i in range(nkeep):
            elistT.append([i,np.real(evalsT[i])])
        elistT.sort(msm_lib.esort)

        # calculate relaxation times from K and T
        tauK = []
        tauT = []
        for i in range(nkeep):
            iiK,lamK = elistK[i]
            if abs(lamK)>1.e-15:
                tauK.append(-1./lamK)
            else:
                tauK.append(0.)
            iT, lamT = elistT[i]
            if lamT>1.e-15:
                tauT.append(-lagt/np.log(lamT))
            else:
                tauT.append(0.)

        # equilibrium probabilities
        ieqK, eK = elistK[0]
        peqK_sum = reduce(lambda x, y: x + y, map(lambda x: rvecsK[x,ieqK],
            range(nkeep)))
        peqK = rvecsK[:,ieqK]/peqK_sum
        ieqT, eT = elistT[0]
        peqT_sum = reduce(lambda x,y: x + y, map(lambda x: rvecsT[x,ieqT],
             range(nkeep)))
        peqT = rvecsT[:,ieqT]/peqT_sum
        return tauT, peqT, rvecsT, lvecsT
#

