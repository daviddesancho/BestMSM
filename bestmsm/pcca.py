import copy
import random
from msm import MSM
import trajectory as traj
import msm_lib
import pcca_lib


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

    def __init__(self, parent, N=2, method="robust"):
        self.parent = parent
        self.N = N
        self.keep_keys = range(N)
        self.macros = self.eigen_group(N=self.N, method=method)

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
        if not hasattr(self.parent, 'lvecsT') and hasattr(self.parent, 'lvecsK'):
            tauT, peqT, self.parent.rvecsT, self.parent.lvecsT = \
                self.parent.tauK, self.parent.peqK, self.parent.rvecsK, self.parent.lvecsK
        elif not hasattr(self.parent, 'lvecsT'):
            tauT, peqT, self.parent.rvecsT, self.parent.lvecsT = \
                   self.parent.calc_eigsT(evecs=True)
        lvecs = self.parent.lvecsT

        # split in desired number of macrostates
        macros = {}
        keep_states = self.parent.keep_states
        keep_keys = self.parent.keep_keys
        macros[0] = range(len(keep_states))
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
        """ Maps trajectory onto the PCCA clusters

        Returns
        -------
        mappedtraj : str
            The mapped trajectory.

        """
        print "\n Mapping trajectory onto macrostates..."
        mappedtraj = []
        keep_states = self.parent.keep_states
        keep_keys = self.parent.keep_keys
        for data in self.parent.data:
            mt = traj.TimeSeries()
            mt.dt = data.dt
            mt.time = data.time
            mt.states = []
            mt.filename = data.filename
            for s in data.states:
                try:
                    mt.states.append([k for k, v in self.macros.iteritems() \
                           if keep_keys.index(s) in v][0])
                except ValueError:
                    #print " not in keep_keys"
                    try:
                        prev = mt.states[-1]
                    except IndexError:
                        pass
            mappedtraj.append(mt)
        self.mappedtraj = mappedtraj
        super(PCCA, self).__init__(self.mappedtraj, keys=range(self.N), lagt=self.parent.lagt)

    def metastability(self):
        """ Calculate metastability according to the definition
        in Chodera et al, J Chem Phys, (2007)

        """
        return pcca_lib.metastability(self.trans)

    def optim(self, nsteps=1, nwrite=None, fout="mc.dat"):
        """ MC optimization using the metastability Q as energy.
        
        Parameters
        ----------
        nsteps : int
            Number of steps per round of MC and per microstate.
        nwrite : int
            Frequency of writing MC output.
        fout : string
            File for output of MC progress.

        Returns
        -------
        macro_opt : dict
            Dictionary with the membership to macrostates.

        """
        print "\n Optimizing the lumped MSM\n"
        out = open(fout, "w")
        out.write("#    iter       q \n")

        nmac = self.N
        nmic = len(self.parent.keep_keys)
        mcsteps = len(self.count)*nsteps*nmic # mc steps per block
        mcsteps_max = nmic*20000 # maximum number of mc steps 
        print self.count
        print self.trans
        q =  self.metastability()
        print " initial:", q
        q_opt = q

        macro = copy.deepcopy(self.macros)
        cont = True
        nmc = 0 # number of mc blocks
        reject = 0
        while cont:
            imc = 0 
            out.write ("%6i %12.10f %10.6e\n"%(imc + nmc*mcsteps,q,1))
            while imc < mcsteps:
                # try ramdom insertion of a microstate in a macrostate
                imac = 0
                jmac = 0
                while imc < mcsteps:
                    imc +=1
                    while True:
                        # choose microstate to move around
                        imic = random.choice(range(nmic))
                        imac = int([x for x in range(nmac) if imic in macro[x]][0])
                        if len(macro[imac]) > 1:
                            # choose destination macrostate
                            jmac = random.choice([x for x in range(nmac) if x is not imac])
                            break
                    # move microstate from i to j
                    macro_new = copy.deepcopy(macro)
                    macro_new[imac].remove(imic)
                    macro_new[jmac].append(imic)
                    # calculate transition count matrix for new mapping
                    count_mac_new = pcca_lib.map_micro2macro(self.parent.count, macro_new, self.parent.keep_states)
                    Tmacro_new = msm_lib.calc_trans(nmac, range(nmac), count_mac_new)
                    # calculate metastability
                    q_new = pcca_lib.metastability(Tmacro_new)
                    delta = pcca_lib.beta(imc,mcsteps)*(q - q_new) # calculate increment (Q is a -Energy)
                    if pcca_lib.metropolis(delta):
                        #print "ACCEPT"
                        macro = copy.deepcopy(macro_new)
                        count_mac = count_mac_new
                        q = q_new
                        if q > q_opt:
                            q_opt = q
                            macro_opt = copy.deepcopy(macro)
                            Tmacro_opt = Tmacro_new
                            self.macro = copy.deepcopy(macro_opt)
                    else:
                        reject+=1
                        #print " REJECT"

                    out.write ("%6i %12.10e %10.6e\n"%(imc + nmc*mcsteps,q,1./pcca_lib.beta(imc,mcsteps)))
                    imc +=1
                cont = False    
        print " final :", q
        print " best :", q_opt
        print " acceptance:",1.-float(reject)/mcsteps

        self.map_trajectory()
        self.do_count()
        self.do_trans()

    def write_mapping(self):
        """ 
        Prints files with the mapping between states and clusters

        """
        for mtraj in self.mappedtraj:
            try:
                idf = mtraj.filename.rfind(".dat")
                filename = mtraj.filename[:idf] + "_mapped_pcca%g.dat"%self.N
            except ValueError:
                filename = mtraj.filename + "_mapped_pcca%g.dat"%self.N
            print " ...writing mapped trajectory at %s"%filename
            fout = open(filename, "w")
            micro_data = [x for x in self.parent.data if x.filename == mtraj.filename][0]
            for x in zip(micro_data.time, micro_data.states, self.data[0].states):
                fout.write("%10.3f %s %8i\n"%(x[0], x[1], x[2]))
            fout.close()
