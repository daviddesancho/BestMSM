#!/usr/bin/env python

import numpy as np
import networkx as nx
import os, math, copy #,numarray,linalg
import itertools,operator
from scipy import linalg as scipyla
import multiprocessing as mp
import cPickle

# thermal energy (kJ/mol)
beta = 1./(8.314e-3*300)

def difference(k1, k2):
    l = len(k1)
    diff = 0
    for i in range(l):
        if k1[i] != k2[i]:
            diff+=1
    return diff

def esort(ei, ej):
    _, eval_i = ei
    _, eval_j = ej

    if eval_j.real > eval_i.real:
        return 1
    elif eval_j.real < eval_i.real:
        return -1
    else:
        return 0

def find_keys(state_keys, trans, manually_remove):
    """ eliminate dead ends """
    keep_states = []
    keep_keys = []
    # eliminate dead ends
    nstate = len(state_keys)
    for i in range(nstate):
        key = state_keys[i]
        summ = 0
        sumx = 0
        for j in range(nstate):
            if j!=i:
                summ += trans[j][i]   # sources
                sumx += trans[i][j] # sinks
        if summ > 0 and sumx > 0 and trans[i][i] > 0 and key not in manually_remove:
            keep_states.append(i)
            keep_keys.append(state_keys[i])
    return keep_states,keep_keys

def connect_groups(keep_states, trans):
    """ check for connected groups """
    connected_groups = []
    leftover = copy.deepcopy(keep_states)
    while len(leftover) > 0:
        #print leftover
        leftover_new = []
        n_old_new_net = 0
        new_net = [ leftover[0] ]
        n_new_net = len(new_net)
        while n_new_net != n_old_new_net:
            for i in range(len(leftover)):
                l = leftover[i]
                if l in new_net:
                    continue
                summ = 0
                for g in new_net:
                    summ += trans[l][g]+trans[g][l]
                if summ > 0:
                    new_net.append(l)
            n_old_new_net = n_new_net
            n_new_net = len(new_net)
            #print " added %i new members" % (n_new_net-n_old_new_net)
        leftover_new = filter(lambda x: x not in new_net, leftover)
        connected_groups.append(new_net)
        leftover = copy.deepcopy(leftover_new)
    return connected_groups

def isnative(native_string, string):
    s = ""
    for i in range(len(string)):
        if string[i]==native_string[i]:
            s+="1"
        else:
            s+="0"
    return s

def mat_mul_v(m, v):
    rows = len(m)
    w = [0]*rows
    irange = range(len(v))
    summ = 0
    for j in range(rows):
        r = m[j]
        for i in irange:
            summ += r[i]*v[i]
        w[j], summ = summ,0
    return w

def dotproduct(v1, v2, sum=sum, imap=itertools.imap, mul=operator.mul):
    return sum(imap(mul,v1,v2))

#def rate_analyze(rate):
#   # calculates eigenvalues and eigenvectors from rate matrix
#   # calculate symmetrized matrix
#   kjisym = kji*(kji.transpose())
#   kjisym = sqrt(kjisym)
#   for j in arange(nstates):
#       kjisym[j,j] = -kjisym[j,j]
#   # calculate eigenvalues and eigenvectors
#   eigvalsym,eigvectsym = linalg.eig(kjisym)
#   # index the solutions
#   index = argsort(-eigvalsym)
#   ieq = index[0]
#   # equilibrium population
#   peq = eigvectsym[:,ieq]**2
#   # order eigenvalues and calculate left and right eigenvectors
#   eigval = zeros((nstates),float)
#   PsiR = zeros((nstates,nstates),float)
#   PsiL = zeros((nstates,nstates),float)       
#   for i in arange(nstates):
#       eigval[i] = eigvalsym[index[i]]
#       PsiR[:,i] = eigvectsym[:,index[i]]*eigvectsym[:,ieq]
#       PsiL[:,i] = eigvectsym[:,index[i]]/eigvectsym[:,ieq]
#   return eigval,PsiR,PsiL,eigvectsym,peq

def propagate(rate, t, pini):
    # propagate dynamics using rate matrix exponential
    expkt = scipyla.expm2(rate*t)
    return mat_mul_v(expkt,pini)

def propagate_eig(elist, rvecs, lvecs, t, pini):
    # propagate dynamics using rate matrix exponential using eigenvalues and eigenvectors 
    nstates = len(pini)
    p = np.zeros((nstates),float)
    for n in range(nstates):
        #print np.exp(-elist[n][1]*t)
        i,e = elist[n]
        p = p + rvecs[:,i]*(np.dot(lvecs[:,i],pini)*\
                np.exp(-abs(e*t)))
    return p

def bootsfiles(traj_list_dt):
    n = len(traj_list_dt)
    traj_list_dt_new = []
    i = 0
    while i < n:
        k = int(np.random.random()*n)
        traj_list_dt_new.append(traj_list_dt[k])
        i += 1
    return traj_list_dt_new 

def boots_pick(filename, blocksize):
    raw = open(filename).readlines()
    lraw = len(raw)
    nblocks = int(lraw/blocksize)
    lblock = int(lraw/nblocks)
    try:
        ib = np.random.randint(nblocks-1)
    except ValueError:
        ib = 0
    return raw[ib*lblock:(ib+1)*lblock]

def onrate(states, target, K, peq):
    # steady state rate
    kon = 0.
    for i in states:
        if i != target:
            if K[target,i] > 0:
                kon += K[target,i]*peq[i]
    return kon

def run_commit(states, K, peq, FF, UU):
    """ calculate committors and reactive flux """
    nstates = len(states)
    # define end-states
    UUFF = UU + FF
    print "   definitely FF and UU states", UUFF
    I = filter(lambda x: x not in UU+FF, states)
    NI = len(I)

    # calculate committors
    b = np.zeros([NI], float)
    A = np.zeros([NI,NI], float)
    for j_ind in range(NI):
        j = I[j_ind]
        summ = 0.
        for i in FF:
            summ += K[i][j]
        b[j_ind] = -summ
        for i_ind in range(NI):
            i = I[i_ind]
            A[j_ind][i_ind] = K[i][j]       
    # solve Ax=b
    Ainv = np.linalg.inv(A)
    x = np.dot(Ainv,b)
    #XX = np.dot(Ainv,A)

    pfold = np.zeros(nstates,float)
    for i in range(nstates):
        if i in UU:
            pfold[i] = 0.0
        elif i in FF:
            pfold[i] = 1.0
        else:
            ii = I.index(i)
            pfold[i] = x[ii]
                        
    # stationary distribution
    pss = np.zeros(nstates,float)
    for i in range(nstates):
        pss[i] = (1-pfold[i])*peq[i]

    # flux matrix and reactive flux
    J = np.zeros([nstates,nstates],float)
    for i in range(nstates):
        for j in range(nstates):
            J[j][i] = K[j][i]*peq[i]*(pfold[j]-pfold[i])

    # dividing line is committor = 0.5 
    #sum_flux = 0
    #left = [x for x in range(nstates) if pfold[x] < 0.5]
    #right = [x for x in range(nstates) if pfold[x] > 0.5]
    #for i in left:
    #    for j in right:
    #        #print "%i --> %i: %10.4e"%(i, j, J[j][i])
    #        sum_flux += J[j][i]

    # dividing line is reaching end states 
    sum_flux = 0
    for i in range(nstates):
        for j in range(nstates):
            if j in FF: #  dividing line corresponds to I to F transitions
                sum_flux += J[j][i]
    print "   reactive flux: %g"%sum_flux

    #sum of populations for all reactant states
    pU = np.sum([peq[x] for x in range(nstates) if pfold[x] < 0.5])
 #   pU = np.sum(peq[filter(lambda x: x in UU, range(nstates))])
    kf = sum_flux/pU
#    print "   binding rate: %g"%kf
    return J, pfold, sum_flux, kf

def calc_count_worker(x):
    states = x[0]
    dt = x[1]
    keys = x[2]
    nkeys = len(keys)
    lagt = x[3]
    sliding = x[4]
    nstates = len(states)
    lag = int(lagt/dt) # number of frames for lag time
    if sliding:
        slider = 1 # every state is initial state
    else:
        slider = lag

    count = np.zeros([nkeys,nkeys], np.int32)
    for i in range(0, nstates-lag, slider):
        j = i + lag
        state_i = states[i]
        state_j = states[j]
        if state_i in keys:
            idx_i = keys.index(state_i)
        if state_j in keys:
            idx_j = keys.index(state_j)
        try:
            count[idx_j][idx_i] += 1
        except UnboundLocalError:
            pass
    return count

def do_boots_worker(x):
    """ Worker function for parallel bootstrapping.

    Parameters:
    ----------
    x : list
        A list containing the trajectory filename, the states, the lag time
        and the total number of transitions.
 
    """

    #print "# Process %s running on input %s"%(mp.current_process(), x[0])
    filetmp, keys, lagt, ncount, slider = x
    nkeys = len(keys)
    finp = open(filetmp, 'rb')
    trans = cPickle.load(finp)
    finp.close()
    ltrans = len(trans)
    # For multiple processes to be independent we seed with pid
    #print 'process id:', os.getpid()
    pid = os.getpid()
    np.random.seed(pid)
    ncount_boots = 0
    count = np.zeros([nkeys, nkeys], np.int32)
    while ncount_boots < ncount:
        itrans = np.random.randint(ltrans)
        count_inp = [trans[itrans][0], trans[itrans][1], keys, lagt, slider]
        c = calc_count_worker(count_inp)
        count += np.matrix(c)
        ncount_boots += np.sum(c)
    D = nx.DiGraph(count)
    #keep_states = sorted(nx.strongly_connected_components(D)[0])
    keep_states = sorted(list(nx.strongly_connected_components(D)), 
                key = len, reverse=True)[0]
    keep_keys = map(lambda x: keys[x], keep_states)
    nkeep = len(keep_keys)
    trans = np.zeros([nkeep, nkeep], float)
    for i in range(nkeep):
        ni = reduce(lambda x, y: x + y, map(lambda x: 
            count[keep_states[x]][keep_states[i]], range(nkeep)))
        for j in range(nkeep):
            trans[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)

    evalsT, rvecsT = scipyla.eig(trans, left=False)
    elistT = []
    for i in range(nkeep):
        elistT.append([i,np.real(evalsT[i])])
    elistT.sort(esort)
    
    tauT = []
    for i in range(1,nkeep):
        _, lamT = elistT[i]
        tauT.append(-lagt/np.log(lamT))
    ieqT, _ = elistT[0]
    peqT_sum = reduce(lambda x,y: x + y, map(lambda x: rvecsT[x,ieqT],
             range(nkeep)))
    peqT = rvecsT[:,ieqT]/peqT_sum
    return tauT, peqT, keep_keys 

def calc_trans(nkeep=None, keep_states=None, count=None):
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
    trans = np.zeros([nkeep, nkeep], float)
    for i in range(nkeep):
        ni = reduce(lambda x, y: x + y, map(lambda x: 
            count[keep_states[x]][keep_states[i]], range(nkeep)))
        for j in range(nkeep):
            trans[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)
    return trans

def calc_rate(nkeep, trans, lagt):
    """ Calculate rate matrix from transition matrix.
    We use the Taylor expansion as described in De Sancho,
    Mittal and Best, J. Chem. Theory Comput. (2013).

    Parameters
    ----------
    nkeep : int
        Number of states in transition matrix.
    trans: np.array
        Transition matrix.
    lagt : float
        The lag time.      

    Returns
    -------
    rate : np.array
        The rate matrix.

    """
    rate = trans/lagt
    
    for i in range(nkeep):
        rate[i][i] = -(np.sum(rate[:i,i]) + np.sum(rate[i+1:,i]))
    return rate

def partial_rate(K, elem):
    """ calculate derivative of rate matrix """
    nstates = len(K[0])
    d_K = np.zeros((nstates,nstates), float)
    for i in range(nstates):
        if i != elem:
            d_K[i,elem] = beta/2.*K[i,elem];
            d_K[elem,i] = -beta/2.*K[elem,i];
    for i in range(nstates):
        d_K[i,i] = -np.sum(d_K[:,i])
    return d_K

def partial_peq(peq, elem):
    """ calculate derivative of equilibrium distribution """
    nstates = len(peq)
    d_peq = []
    for i in range(nstates):
        if i != elem:
            d_peq.append(beta*peq[i]*peq[elem])
        else:
            d_peq.append(-beta*peq[i]*(1.-peq[i]))
    return d_peq

def partial_pfold(states, K, d_K, FF, UU, elem):
    """ calculate derivative of pfold """
    nstates = len(states)
    # define end-states
    UUFF = UU+FF
    I = filter(lambda x: x not in UU+FF, range(nstates))
    NI = len(I)
    # calculate committors
    b = np.zeros([NI], float)
    A = np.zeros([NI,NI], float)
    db = np.zeros([NI], float)
    dA = np.zeros([NI,NI], float)
    for j_ind in range(NI):
        j = I[j_ind]
        summ = 0.
        sumd = 0.
        for i in FF:
            summ += K[i][j]
            sumd+= d_K[i][j]
        b[j_ind] = -summ
        db[j_ind] = -sumd
        for i_ind in range(NI):
            i = I[i_ind]
            A[j_ind][i_ind] = K[i][j]
            dA[j_ind][i_ind] = d_K[i][j]

    # solve Ax + Bd(x) = c
    Ainv = np.linalg.inv(A)
    pfold = np.dot(Ainv,b)
    x = np.dot(Ainv,db - np.dot(dA,pfold))

    dpfold = np.zeros(nstates,float)
    for i in range(nstates):
        if i in UU:
            dpfold[i] = 0.0
        elif i in FF:
            dpfold[i] = 0.0
        else:
            ii = I.index(i)
            dpfold[i] = x[ii]
    return dpfold

def partial_flux(states,peq,K,pfold,d_peq,d_K,d_pfold,target):
    """ calculate derivative of flux """
    # flux matrix and reactive flux
    nstates = len(states)
    sum_d_flux = 0
    d_J = np.zeros((nstates,nstates),float)
    for i in range(nstates):
        for j in range(nstates):
            d_J[j][i] = d_K[j][i]*peq[i]*(pfold[j]-pfold[i]) + \
                K[j][i]*d_peq[i]*(pfold[j]-pfold[i]) + \
                K[j][i]*peq[i]*(d_pfold[j]-d_pfold[i])
            if j in target and K[j][i]>0: #  dividing line corresponds to I to F transitions                        
                sum_d_flux += d_J[j][i]
    return sum_d_flux

def propagate_worker(x):
    """ propagate dynamics using rate matrix exponential"""
    rate, t, pini = x
    expkt = scipyla.expm2(rate*t)
    popul = mat_mul_v(expkt, pini)
    return popul 

def gen_path_lengths(keys, J, pfold, flux, FF, UU):
    """ use BHS prescription for defining path lenghts """
    nkeys = len(keys)
    I = [x for x in range(nkeys) if x not in FF+UU]
    Jnode = []
    # calculate flux going through nodes
    for i in range(nkeys):
        Jnode.append(np.sum([J[i,x] for x in range(nkeys) \
                             if pfold[x] < pfold[i]]))
    # define matrix with edge lengths
    Jpath = np.zeros((nkeys, nkeys), float)
    for i in UU:
        for j in I + FF:
            if J[j,i] > 0:
                Jpath[j,i] = np.log(flux/J[j,i]) + 1
    for i in I:
        for j in [x for x in FF+I if pfold[x] > pfold[i]]:
            if J[j,i] > 0:
                Jpath[j,i] = np.log(Jnode[j]/J[j,i]) + 1e-99
    return Jnode, Jpath
