""" 
================
 MSM 
================
"""

import os
import sys
import tempfile
import cPickle
import numpy as np
import networkx as nx
import scipy.linalg as scipyla
import matplotlib.pyplot as plt
import multiprocessing as mp
import msm_lib
import visual_lib

class MasterMSM(object):
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

        Returns: 
        -------
        msm : dict
            A dictionary with the multiple instances of the MSM class.

        """

        # defining lag times to produce the MSM
        lagtimes = self.dt*np.array(range(1,20,2))

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
                ax.plot(data[:,0], data[:,n+1], label=n)
            ax.set_xlabel(r'Time', fontsize=16)
            ax.set_ylabel(r'$\tau$', fontsize=16)
            plt.show()

        return msms

#    def do_pcca(self, lagt=10, N=2, optim=True):
#        """ Do PCCA clustering
#
#        Parameters:
#        -----------
#        lagt : float
#            The lag time.
#
#        N : int
#            The number of clusters.
#
#        optim : bool
#            Whether optimization of the clustering is desired.
#
#        """
#
#        return self.msms[lagt].pcca(lagt=lagt, N=N, optim=optim)

class MSM(object):
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

    def __init__(self, data, keys=None, lagt=None, evecs=False):
        # merge state keys from all trajectories
        self.keys = keys
        self.data = data
        self.lagt = lagt
        self.count = self.calc_count_multi(lagt)
        self.keep_states, self.keep_keys = self.check_connect()
        self.trans = self.do_trans(lagt)
        self.rate = self.do_rate(lagt)
        print self.trans
        print self.rate
        if not evecs:
            self.tauT, self.peqT = self.calc_eigs(lagt=lagt)
        else:
            self.tauT, self.peqT, self.rvecsT, self.lvecsT = \
                    self.calc_eigs(lagt=lagt, evecs=True)

    def calc_count_multi(self, lagt=None, nproc=None):
        """ Calculate transition count matrix in parallel 
        
        """

        print "\n Calculating transition count matrix..."
        if not nproc:           
            nproc = mp.cpu_count()
            if len(self.data) < nproc:
                nproc = len(self.data)
                print "\n    ...running on %g processors"%nproc
        pool = mp.Pool(processes=nproc)
        mpinput = map(lambda x: [x.states, x.dt, self.keys, lagt], self.data)
        result = pool.map(msm_lib.calc_count_worker, mpinput)
        pool.close()
        pool.join()

        count = reduce(lambda x, y: np.matrix(x) + np.matrix(y), result)
        return np.array(count)
 
    def calc_count_seq(self, lagt=None):
        """ Calculate transition count matrix sequentially 
        
        """

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

    def do_trans(self, lagt=None):
        """ Wrapper for transition matrix calculation.

        Parameters:
        ----------
        lagt : float
            Lag time for construction of MSM.

        Returns:
        -------
        array
            The transition probability matrix.    
        
        """

        print "\n Calculating transition matrix ..."
        nkeep = len(self.keep_states)
        keep_states = self.keep_states
        count = self.count
        return msm_lib.calc_trans(nkeep, keep_states, count)
    
    def do_rate(self, lagt=None):
        """ Wrapper for rate calculation using the msm_lib.calc_rate 
        function.

        Parameters
        ----------
        lagt : float
            The lag time.

        Returns
        -------
        array
            The rate matrix as calculated by msm_lib.calc_rate

        """
        print "\n Calculating rate matrix ..."
        nkeep = len(self.keep_states)
        return msm_lib.calc_rate(nkeep, self.trans, lagt)

    def calc_eigs(self, lagt=None, evecs=False):
        """ Calculate eigenvalues and eigenvectors
        
        Parameters:
        -----------
        lagt : float
            Lag time used for constructing MSM.
        evecs : bool
            Whether we want eigenvectors or not.

        Returns:
        -------
        tauT : numpy array
            Relaxation times from T.
        peqT : numpy array
            Equilibrium probabilities from T.
        rvecsT : numpy array, optional
            Right eigenvectors of T, sorted.
        
        lvecsT : numpy array, optional
            Left eigenvectors of T, sorted.

        """
        print "\n Calculating eigenvalues and eigenvectors"
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
        for i in range(1, nkeep):
            iiK, lamK = elistK[i]
            tauK.append(-1./lamK)
            iiT, lamT = elistT[i]
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

        if not evecs:
            return tauT, peqT
        else:
            # sort eigenvectors
            rvecsT_sorted = np.zeros((nkeep, nkeep), float)
            lvecsT_sorted = np.zeros((nkeep, nkeep), float)
            for i in range(nkeep):
                iiT, lamT = elistT[i]
                rvecsT_sorted[:,i] = rvecsT[:,iiT]
                lvecsT_sorted[:,i] = lvecsT[:,iiT]
            return tauT, peqT, rvecsT_sorted, lvecsT_sorted

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
            nboots = 100
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
        pool.close()
        pool.join()

        tauT_boots = [x[0] for x in result]
        peqT_boots = [x[1] for x in result]
        keep_keys_boots = [x[2] for x in result]
        # TODO. rewrite so that auxiliary functions are in msm_lib
        tau_ave = []
        tau_std = []
        tau_keep = []
        for n in range(len(self.keys)-1):
            try:
                data = [x[n] for x in tauT_boots if not np.isnan(x[n])]
                tau_ave.append(np.mean(data))
                tau_std.append(np.std(data))
                tau_keep.append(data)
            except IndexError:
                continue
        peq_ave = []
        peq_std = []
        peq_indexes = []
        peq_keep = []
        for k in self.keys:
            peq_indexes.append([x.index(k) if k in x else None for x in keep_keys_boots])

        for k in self.keys:
            l = self.keys.index(k)
            data = []
            for n in range(nboots):
                if peq_indexes[l][n] is not None:
                    data.append(peqT_boots[n][peq_indexes[l][n]])
            try:
                peq_ave.append(np.mean(data))
                peq_std.append(np.std(data))
                peq_keep.append(data)
            except RuntimeWarning:
                peq_ave.append(0.)
                peq_std.append(0.)

        if plot:
            data = np.array(data)
            fig = plt.figure()
            fig.set_facecolor('white')
            ax = fig.add_subplot(2,1,1)
            binning = np.arange(0, 1, 0.01)
            for n in range(10):
                ax.hist(peq_keep[n], bins=binning)
            ax.set_xlabel(r'$P_{eq}$')

            ax = fig.add_subplot(2,1,2)
            binning = np.arange(np.log10(np.min(tau_ave)) - 1,np.log10(np.max(tau_ave)) + 1 , 0.05)
            for n in range(10):
                ax.hist(np.log10(tau_keep[n]), bins=binning)
            ax.set_xlabel(r'$\tau$')
            plt.show()

        return tau_ave, tau_std, peq_ave, peq_std

#    def pcca(self, lagt=None, N=2, optim=False):
#        """ Wrapper for carrying out PCCA clustering
#
#        Parameters
#        ----------
#        lagt : float
#            The lag time.
#
#        N : int
#            The number of clusters.
#
#        optim : bool
#            Whether optimization is desired.
#        
#        """
#
#        return pcca.PCCA(parent=self, lagt=lagt, optim=optim)


    def do_pfold(self, FF=None, UU=None, dot=False):
        """ Wrapper to calculate reactive fluxes and committors using the 
        Berzhkovskii-Hummer-Szabo method, J Chem Phys (2009)
        
        Parameters
        ----------
        lagt : float
            The lag time.
        FF : list
            Folded states.
        UU : list
            Unfolded states.
        dot : string
            Filename to output dot graph.

        Returns
        -------
        J : array
            The flux matrix.
        pfold : list
            The values of the committor.
        kf : float
            The foldign rate UU -> FF

        """
        print "\n Calculating commitment probabilities and fluxes..."
        _states = range(len(self.keep_states))
        if isinstance(FF, list):
            _FF = [self.keep_states.index(x) for x in FF]
        else:
            _FF = [self.keep_states.index(FF)]
        if isinstance(UU, list):
            _UU = [self.keep_states.index(x) for x in UU]
        else:
            _UU = [self.keep_states.index(UU)]

        # do the calculation
        J, pfold, kf = msm_lib.run_commit(_states, self.rate, self.peqT, _FF, _UU)

        # write graph in dot format
        if dot:
            D = nx.DiGraph(J)
#            visual_lib.write_dot(D, nodeweight=self.peqT, edgeweight=self.trans, out="out.dot")
            visual_lib.write_dot(D, nodeweight=self.peqT, out="out.dot")
#            visual_lib.write_dot(D, out="out.dot")

        return J, pfold, kf
