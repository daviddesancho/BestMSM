""" 
================
 MSM 
================
"""

import os
import numpy as np
import scipy.linalg as scipyla
import msm_lib
import tempfile
import cPickle
import networkx as nx
import matplotlib.pyplot as plt
import multiprocessing as mp

class MasterMSM:
    """
    A class for constructing the MSM

    Parameters
    ----------
    trajectories : list of str
        Set of trajectories used for the construction.

    filekeys : str
        A file from which the states are read. If not defined, then
        keys are automatically generated from the time series.


    Attributes
    ----------
    keys : list of str
        Names of states.

    count : np array
        Transition count matrix.

    keep_states : list of str
        Names of states after removing not strongly connected sets.

    msms : dict
        A dictionary containing MSMs for different lag times.
        
    """

    def __init__(self, trajectories, filekeys=None):
        self.data = trajectories
        try:
            self.keys = map(lambda x: x.split()[0],
                    open(filekeys, "r").readlines())
        except TypeError:
            self.keys = self.__merge_trajs()
        self.dt = self.__min_dt()
        self.__out()
        self.msms = {}

    def __merge_trajs(self):
        """ Merge all trajectories into a consistent set.
        
        """
        new_keys = []
        for traj in self.data:
            new_keys += filter(lambda x: x not in new_keys, traj.keys)
        return new_keys

    def __min_dt(self):
        """ Find minimum dt in trajectories.
        
        """
        return np.min(map(lambda x: x.dt, self.data))

    def __out(self):
        """ Output description to user.
        
        """
        print "\n Building MSM from \n",map(lambda x: x.filename, self.data)
        print "     # states: %g"%(len(self.keys))
        print self.keys

    def do_msm(self, lagt):
        """ Construct MSM for specific value of lag time.
        
        Parameters:
        -----------
        lagt : float
            The lag time.

        """

        self.msms[lagt] = MSM(self.data, self.keys, lagt)


    def chapman_kolmogorov(self, plot=True, N=1):
        """ Carry out Chapman-Kolmogorov test.

        Parameters:
        -----------
        plot : bool
            Whether the lag time dependence of the relaxation times should be plotted.

        N: int
            The number of modes for which we are building the MSM.

        Returns : dict
            A dictionary with the multiple instances of the MSM class.

        """

        # defining lag times to produce the MSM
        lagtimes = self.dt*np.array([1, 2, 5, 10, 20, 50, 100, 200, 500])

        # create MSMs at multiple lag times
        msms = {}
        data = []
        for lagt in lagtimes:
            print "\n Generating MSM at lag time: %g"%lagt
            msms[lagt] = MSM(self.data, self.keys, lagt)
            print "\n    Count matrix:\n", msms[lagt].count
            print "\n    Transition matrix:\n", msms[lagt].trans
            dat = [lagt]
            for n in range(N):
                dat.append(msms[lagt].tauT[n])
            data.append(dat)

        if plot:
            data = np.array(data)
            fig, ax = plt.subplots(facecolor='white')
            for n in range(N):
                ax.loglog(data[:,0], data[:,n+1], label=n)
            ax.set_xlabel(r'Time', fontsize=16)
            ax.set_ylabel(r'$\tau$', fontsize=16)
            plt.show()

        return msms

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
        self.lagt = lagt
        self.count = self.calc_count_multi(lagt)
        self.keep_states, self.keep_keys = self.check_connect()
        self.trans = self.calc_trans(lagt)
        self.rate = self.calc_rate(lagt)
        self.tauT, self.peqT, self.rvecsT, self.lvecsT = \
            self.calc_eigs(lagt)

    def calc_count_multi(self, lagt=None, nproc=None):
        """ Calculate transition count matrix in parallel 
        
        """
        if not nproc:           
            nproc = mp.cpu_count()
            if len(self.data) < nproc:
                nproc = len(self.data)
                print "\n    ...running on %g processors"%nproc
        pool = mp.Pool(processes=nproc)
        mpinput = map(lambda x: [x.states, x.dt, self.keys, lagt], self.data)
        result = pool.map(msm_lib.calc_count_worker, mpinput)
        count = reduce(lambda x, y: np.matrix(x) + np.matrix(y), result)
        return np.array(count)
 
    def calc_count_seq(self, lagt=None):
        """ Calculate transition count matrix sequentially 
        
        """
        print self.keys
        keys = self.keys
        nkeys = len(keys)
        count = np.zeros([nkeys,nkeys], int)

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
        """ Check connectivity of rate matrix. 
        
        """
        print "\n    ...checking connectivity:"
        D = nx.DiGraph(self.count)
        keep_states = sorted(nx.strongly_connected_components(D)[0])
        keep_keys = map(lambda x: self.keys[x], keep_states)
        print "          %g states in largest subgraph"%len(keep_keys)
        return keep_states, keep_keys

    def calc_trans(self, lagt=None):
        """ Calculate transition matrix.

        Parameters:
        ----------
        lagt : float
            Lag time for construction of MSM.

        Returns:
        -------
        trans : array
            The transition probability matrix.    
        
        """
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
        """ Calculate rate matrix using a Taylor series as
        described in De Sancho et al J Chem Theor Comput (2013)
        :param lagt:
        """
        nkeep = len(self.keep_states)
        rate = self.trans/lagt
        for i in range(nkeep):
            rate[i][i] = -(np.sum(rate[:i,i]) + np.sum(rate[i+1:,i]))
        return rate

    def calc_eigs(self, lagt=None):
        """ Calculate eigenvalues and eigenvectors
        :param lagt:
        """

        evalsK, lvecsK, rvecsK = \
                   scipyla.eig(self.rate, left=True)
        evalsT, lvecsT, rvecsT = \
                scipyla.eig(self.trans, left=True)
#        # sort modes
        nkeep = len(self.keep_states)
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
        for i in range(1,nkeep):
            iiK,lamK = elistK[i]
            tauK.append(-1./lamK)
            iT, lamT = elistT[i]
            tauT.append(-lagt/np.log(lamT))

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

    def boots(self, nboots=None, nproc=None, plot=False):
        """ Bootstrap the simulation data to calculate errors

        Parameters:
        -----------
        nboots : int
            Number of bootstrap samples

        nproc : int
            Number of processors to use

        plot : bool
            Whether we want to plot the distribution of tau / peq
        
        Returns:
        --------
        tau_err : array
            Errors for the relaxation times

        peq_err : array
            Errors for the equilibrium probabilities

        """

        print "\n Doing bootstrap tests:"
        # how much data is here?
        # generate trajectory list for easy handling 
        trajs = [[x.states, x.dt] for x in self.data]
        ltraj = [len(x[0])*x[1] for x in trajs]
        timetot = np.sum(ltraj) # total simulation time
        ncount = np.sum(self.count)
        print "     Total time: %g"%timetot
        print "     Number of trajectories: %g"%len(trajs)
        print "     Total number of transitions: %g"%ncount

        # how many resamples?
        if not nboots:
            nboots = 10
        print "     Number of resamples: %g"%nboots

        # how many trajectory fragments?
        ltraj_median = np.median(ltraj)
        if ltraj_median > timetot/10.:
            print "     ...splitting trajectories"
        while ltraj_median > timetot/10.:
            trajs_new = []
            #cut trajectories in chunks
            for x in trajs:
                lx = len(x[0])
                trajs_new.append([x[0][:lx/2], x[1]])
                trajs_new.append([x[0][lx/2:], x[1]])
            trajs = trajs_new
            ltraj = [len(x[0])*x[1] for x in trajs]
            ltraj_median = np.median(ltraj)

        # save trajs
        fd, filetmp = tempfile.mkstemp()
        file = os.fdopen(fd, 'wb')   
        cPickle.dump(trajs, file, protocol=cPickle.HIGHEST_PROTOCOL)
        file.close()

        print "     Number of trajectories: %g"%len(trajs)
        print "     Median of trajectory length: %g"%ltraj_median

        # now do it
        print "     ...doing bootstrap analysis"
        # multiprocessing options
        if not nproc:           
            nproc = mp.cpu_count()
        print "     ...running on %g processors"%nproc
        pool = mp.Pool(processes=nproc)
        multi_boots_input = map(lambda x: [filetmp, self.keys, self.lagt, ncount], 
                range(nboots))
        # TODO: find more elegant way to pass arguments
        result = pool.map(msm_lib.do_boots_worker, multi_boots_input)
        tauT_boots = [x[0] for x in result]
        peqT_boots = [x[1] for x in result]
        keep_keys_boots = [x[2] for x in result]
        # TODO. rewrite so that auxiliary functions are in msm_lib
        tau_err = []
        for n in range(len(self.keys)-1):
            tau_err.append(np.std([x[n] for x in tauT]))
        peq_err = []
        for n in range(len(self.keys)):
            print zip(peqT, keep_keys_boots)
