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
        if not hasattr(self.parent, 'lvecsT'):
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
        print self.macros
        keep_states = self.parent.keep_states
        keep_keys = self.parent.keep_keys
        for data in self.parent.data:
            mt = traj.TimeSeries()
            mt.dt = data.dt
            mt.time = data.time
            mt.states = []
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

    def optim(self, fout="mc.dat"):
        """ MC optimization using the metastability Q as energy.
        
        Parameters
        ----------
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
        nmic = len(self.parent.keys)
        mcsteps = len(self.count)*1000*nmic # mc steps per block
        mcsteps_max = nmic*20000 # maximum number of mc steps 
        q =  self.metastability()
        print " initial:", q
        q_opt = q

        macro = copy.deepcopy(self.macros)
        cont = True
        nmc = 0 # number of mc blocks
        reject = 0
        while cont:
            imc = 0 
            out.write ("%6i %12.10f\n"%(imc + nmc*mcsteps,q))
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
                    count_mac_new = pcca_lib.map_micro2macro(self.parent.count, macro_new)
                    Kmacro_new,Tmacro_new = msm_lib.calc_rate(nmac, count_mac_new, self.lagt)
#                    # calculate metastability
#                    q_new = metastability(Tmacro_new)
#                    #print "Q new: %g"%q_new
#                    #print "temp: %g ="%temp
#                    #print " imc= %g; beta = %g"%(imc,beta)
#                    delta = beta(imc,mcsteps)*(q - q_new) # calculate increment (Q is a -Energy)
#                    #print delta
#                    if metropolis(delta):
#                        #print "ACCEPT"
#                        macro = copy.deepcopy(macro_new)
#                        count_mac = count_mac_new
#                        q = q_new
#                        if q > q_opt:
#                            q_opt = q
#                            macro_opt = copy.deepcopy(macro)
#                    else:
#                        reject+=1
#                        #print " REJECT"
#
#                    if imc%100==0:
#                        out.write ("%6i %12.10e %12.10e\n"%(imc + nmc*mcsteps,q,1./beta(imc,mcsteps)))
#                nmc +=1
#            cont = False    
#        print " final :",q_opt
#        print macro_opt
#        print " acceptance:",1.-float(reject)/mcsteps
#        return macro_opt
#            return pcca_lib.do_mc(self.macros, count=self.parent.count, lagt=self.parent.lagt)
