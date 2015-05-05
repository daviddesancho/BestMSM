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

    def do_msm(self, lagt, sliding=True):
        """ Construct MSM for specific value of lag time.
        
        Parameters:
        -----------
        lagt : float
            The lag time.

        """
        self.msms[lagt] = MSM(self.data, keys=self.keys, lagt=lagt)
        self.msms[lagt].do_count(sliding=sliding)
        self.msms[lagt].do_trans()

    def chapman_kolmogorov(self, plot=True, N=1, sliding=True, error=True):
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
        lagtimes = self.dt*np.array([1] + range(5,50,5))

        # create MSMs at multiple lag times
        self.msms = {}
        for lagt in lagtimes:
            print "\n Generating MSM at lag time: %g"%lagt
            self.msms[lagt] = MSM(self.data, self.keys, lagt, sliding=sliding)
            print "\n    Count matrix:\n", self.msms[lagt].count
            print "\n    Transition matrix:\n", self.msms[lagt].trans
            if error:               
                tau_ave, tau_std, peq_ave, peq_std = self.msms[lagt].boots(nboots=48)
                self.msms[lagt].tau_std = tau_std
                self.msms[lagt].tau_ave = tau_ave
                self.msms[lagt].peq_std = peq_std
                self.msms[lagt].peq_ave = peq_std

        if plot:
            fig, ax = plt.subplots(facecolor='white')
            for n in range(N):
                data = [self.msms[x].tauT[n] for x in lagtimes]
                if not error:
                    ax.plot(lagtimes, data, 'o-', label=n)
                else:
                    ebar = [self.msms[x].tau_std[n] for x in lagtimes]
                    ax.errorbar(lagtimes, data, yerr=ebar, fmt='o', label=n)
            ax.set_xlabel(r'Time', fontsize=16)
            ax.set_ylabel(r'$\tau$', fontsize=16)
            plt.show()

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
    data : list of str
        Set of trajectories used for the construction.

    keys : list of str
        Set of states in the model.

    lag : float
        Lag time for building the MSM.

    sliding : bool
        Whether the method of sliding windows or independent counts 
        are used.

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
        print keys

    def do_count(self, sliding=True):
        self.count = self.calc_count_multi(sliding=sliding)
        self.keep_states, self.keep_keys = self.check_connect()

    def do_trans(self, evecs=False):
        """ 
        Wrapper for transition matrix calculation.
        
        """
        print "\n Calculating transition matrix ..."
        nkeep = len(self.keep_states)
        keep_states = self.keep_states
        count = self.count
        self.trans = msm_lib.calc_trans(nkeep, keep_states, count)
        if not evecs:
            self.tauT, self.peqT = self.calc_eigsT()
        else:
            self.tauT, self.peqT, self.rvecsT, self.lvecsT = \
                    self.calc_eigsT(evecs=True)

    def do_rate(self, evecs=False):
        """ Wrapper for rate calculation using the msm_lib.calc_rate 
        function.

        """
        print "\n Calculating rate matrix ..."
        nkeep = len(self.keep_states)
        self.rate = msm_lib.calc_rate(nkeep, self.trans, self.lagt)
        print self.rate
        if not evecs:
            self.tauK, self.peqK = self.calc_eigsK()
        else:
            self.tauK, self.peqK, self.rvecsK, self.lvecsK = \
                    self.calc_eigsK(evecs=True)

    def calc_count_multi(self, sliding=True, nproc=None):
        """ Calculate transition count matrix in parallel 
        
        """
        print "\n Calculating transition count matrix...\n"
        print "   sliding window option: ", sliding

        # define multiprocessing options
        if not nproc:           
            nproc = mp.cpu_count()
            if len(self.data) < nproc:
                nproc = len(self.data)
                print "\n    ...running on %g processors"%nproc
        pool = mp.Pool(processes=nproc)
        # generate multiprocessing input
        mpinput = map(lambda x: [x.states, x.dt, self.keys, self.lagt, sliding], self.data)
        # run counting using multiprocessing
        result = pool.map(msm_lib.calc_count_worker, mpinput)
        pool.close()
        pool.join()
        # add up all independent counts
        count = reduce(lambda x, y: np.matrix(x) + np.matrix(y), result)
        return np.array(count)
 
    def calc_count_seq(self, sliding=True):
        """ Calculate transition count matrix sequentially 
        
        """
        keys = self.keys
        nkeys = len(keys)
        count = np.zeros([nkeys,nkeys], int)

        for traj in self.data:
            raw = traj.states
            dt = traj.dt
            if sliding:
                lag = 1 # every state is initial state
            else:
                lag = int(self.lagt/dt) # number of frames for lag time
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
        print self.count
        D = nx.DiGraph(self.count)
        keep_states = sorted(list(nx.strongly_connected_components(D)), 
                key = len, reverse=True)[0]
        keep_states.sort()
        #keep_states = sorted(nx.strongly_connected_components(D)[0])
        keep_keys = map(lambda x: self.keys[x], keep_states)
        print "          %g states in largest subgraph"%len(keep_keys)
        return keep_states, keep_keys

    def calc_eigsK(self, evecs=False):
        """ Calculate eigenvalues and eigenvectors of rate matrix K
        
        Parameters:
        -----------
        evecs : bool
            Whether we want eigenvectors or not.

        Returns:
        -------
        tauK : numpy array
            Relaxation times from K.
        peqK : numpy array
            Equilibrium probabilities from K.
        rvecsK : numpy array, optional
            Right eigenvectors of K, sorted.
        lvecsK : numpy array, optional
            Left eigenvectors of K, sorted.

        """
        print "\n Calculating eigenvalues and eigenvectors of K"
        evalsK, lvecsK, rvecsK = \
                   scipyla.eig(self.rate, left=True)
        # sort modes
        nkeep = len(self.keep_states)
        elistK = []
        for i in range(nkeep):
            elistK.append([i,np.real(evalsK[i])])
        elistK.sort(msm_lib.esort)

        # calculate relaxation times from K and T
        tauK = []
        for i in range(1, nkeep):
            iiK, lamK = elistK[i]
            tauK.append(-1./lamK)

        # equilibrium probabilities
        ieqK, eK = elistK[0]
        peqK_sum = reduce(lambda x, y: x + y, map(lambda x: rvecsK[x,ieqK],
            range(nkeep)))
        peqK = rvecsK[:,ieqK]/peqK_sum
        if not evecs:
            return tauK, peqK
        else:
            # sort eigenvectors
            rvecsK_sorted = np.zeros((nkeep, nkeep), float)
            lvecsK_sorted = np.zeros((nkeep, nkeep), float)
            for i in range(nkeep):
                iiK, lamK = elistK[i]
                rvecsK_sorted[:,i] = rvecsK[:,iiK]
                lvecsK_sorted[:,i] = lvecsK[:,iiK]
            return tauK, peqK, rvecsK_sorted, lvecsK_sorted

    def calc_eigsT(self, evecs=False):
        """ Calculate eigenvalues and eigenvectors of transition matrix T
        
        Parameters:
        -----------
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
        print "\n Calculating eigenvalues and eigenvectors of T"
        evalsT, lvecsT, rvecsT = \
                scipyla.eig(self.trans, left=True)
        # sort modes
        nkeep = len(self.keep_states)
        elistT = []
        for i in range(nkeep):
            elistT.append([i,np.real(evalsT[i])])
        elistT.sort(msm_lib.esort)

        # calculate relaxation times 
        tauT = []
        for i in range(1, nkeep):
            iiT, lamT = elistT[i]
            tauT.append(-self.lagt/np.log(lamT))

        # equilibrium probabilities
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

    def boots(self, nboots=None, nproc=None, plot=False, slider=False):
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

        multi_boots_input = map(lambda x: [filetmp, self.keys, self.lagt, ncount, 
            slider], range(nboots))
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

    def rate_scale(self, iscale=None, scaling=1):
        """ Scaling columns of the rate matrix by a certain value
        
        Parameters
        ----------
        states : list
            The list of state keys we want to scale.
        val : float
            The scaling factor.

        """
        print "\n scaling kiS by %g"%scaling
        nkeep = len(self.keep_keys)
        jscale = [self.keep_keys.index(x) for x in iscale]
        for j in jscale:
            for l in range(nkeep):
                self.rate[l,j] = self.rate[l,j]*scaling
#        tauK, peqK, rvecsK_sorted, lvecsK_sorted = self.rate.calc_eigs(lagt=lagt, evecs=False)


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
            _FF = [self.keep_keys.index(x) for x in FF]
        else:
            _FF = [self.keep_keys.index(FF)]
        if isinstance(UU, list):
            _UU = [self.keep_keys.index(x) for x in UU]
        else:
            _UU = [self.keep_keys.index(UU)]

        # do the calculation
        J, pfold, sum_flux, kf = msm_lib.run_commit(_states, \
                self.rate, self.peqK, _FF, _UU)

        # write graph in dot format
        if dot:
            D = nx.DiGraph(J)
#            visual_lib.write_dot(D, nodeweight=self.peqT, edgeweight=self.trans, out="out.dot")
            visual_lib.write_dot(D, nodeweight=self.peqT, out="out.dot")
#            visual_lib.write_dot(D, out="out.dot")

        return J, pfold, sum_flux, kf

    def sensitivity(self, FF=None, UU=None, dot=False):
        """ Sensitivity analysis of the states in the network.
        De Sancho, Kubas, Blumberger and Best (In preparation, 
        2014)

        Parameters
        ----------
        FF : list
            Folded states.
        UU : list
            Unfolded states.

        Returns
        -------
        dJ : list
            Derivative of flux.
        d_peq : list
            Derivative of equilibrium populations.
        d_kon : list
            Derivative of global rate.
        kon : float
            Global rate.

        """
        nkeep = len(self.keep_states)
        K = self.rate
        peqK = self.peqK

        # calculate pfold 
        J, pfold, sum_flux, kon = self.do_pfold(FF=FF, UU=UU)

        # 
        if isinstance(FF, list):
            _FF = [self.keep_keys.index(x) for x in FF]
        else:
            _FF = [self.keep_keys.index(FF)]
        if isinstance(UU, list):
            _UU = [self.keep_keys.index(x) for x in UU]
        else:
            _UU = [self.keep_keys.index(UU)]

        pu = np.sum([peqK[x] for x in range(nkeep) if x not in _FF])
        dJ = []
        d_pu = []
        d_kon = []
        for s in range(nkeep):
            d_K = msm_lib.partial_rate(K, s)
            d_peq = msm_lib.partial_peq(peqK, s)
            d_pfold = msm_lib.partial_pfold(range(nkeep), K, d_K, _FF, _UU, s)
            dJ.append(msm_lib.partial_flux(range(nkeep), peqK, K, pfold, \
                    d_peq, d_K, d_pfold, _FF))
            d_pu.append(np.sum([d_peq[x] for x in range(nkeep) \
                    if x not in _FF]))
            d_kon.append((dJ[-1]*pu - sum_flux*d_pu[-1])/pu**2)

        return dJ, d_peq, d_kon, kon 

    def propagateK(self, p0=None, init=None):
        """ Propagation of rate matrix using matrix exponential 
        
        Parameters
        ----------
        p0 : string
            Filename with initial population.
        init : string
            State name corresponding to 100% initial population.
        
        Returns
        -------
        popul : array
            Population of all states as a function of time.
        pnorm : array
            Population of all states as a function of time - normalized.
        """
         # time for relaxation
        logt = np.arange(np.log10(self.lagt), 7., 0.25)
        time = 10**logt
        ltime = len(time)

        nkeep = len(self.keep_states)
        if p0 is not None:
            try:
                print " reading initial population from file: %s"%p0
                pini = [float(y) for y in \
                        filter(lambda x: x.split()[0] not in ["#","@"],
                        open(p0, "r").readlines())]
            except TypeError:
                print " p0 is not file"
                print " exiting here"
                return
        elif init is not None:
            print " initializing all population in states"
            print init
            pini = [self.peqK[x] if self.keep_keys[x] in init else 0. for x in range(nkeep)]
        # check normalization and size
        if len(pini) != nkeep:
            print " initial population vector and state space have different sizes"
            print " stopping here" 
            return
        else:
            sum_pini = np.sum(pini)
            pini_norm = [np.float(x)/sum_pini for x in pini]

        # propagate rate matrix : parallel version
        nproc = mp.cpu_count()
        pool = mp.Pool(processes=nproc)
        pool_input = [(self.rate, t, pini_norm) for t in time]
        popul = pool.map(msm_lib.propagate_worker, tuple(pool_input))
        pool.close()
        pool.join()

        ## normalize relaxation
        #imax = np.argmax(ptot)
        #maxpop = ptot[imax]
        #imin = np.argmin(ptot)
        #minpop = ptot[imin]
        #if imax < imin:
        #    pnorm = map(lambda x: (x-minpop)/(maxpop-minpop), ptot)
        #else:
        #    pnorm = map(lambda x: 1 - (x-minpop)/(maxpop-minpop), ptot)
        return time, popul #popul #, pnorm
