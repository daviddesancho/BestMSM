import sys,copy,itertools,random,math
import numpy as np
from scipy import linalg as scipyla

"""useful functions for clustering"""

def calc_trans(mat):
    # calculates transition matrix given count matrix
    n = len(mat)
    trans = np.zeros((n,n),float)
    for i in range(n):
        ni = reduce(lambda x,y:x+y, map(lambda x: mat[x,i], range(n)))
        for j in range(n):
            trans[j,i] = float(mat[j,i])/float(ni)
    return trans

def calc_rate(mat,t):
    # calculate rate matrix
    n = len(mat)
    T = calc_trans(mat)
    K = T/t
    for i in range(n):
        K[i][i] = -(np.sum(K[:i,i]) + np.sum(K[i+1:,i]))
    return K,T

def map_micro2macro(cmic,mac):
    # maps microstates into macrostates and returns count matrix
    n = len(cmic)
    m = len(mac)
    cmac = np.zeros((m,m),int)
    for i in range(m):
        for j in range(m):
            cmac[j,i] = reduce(lambda x,y: x+y,map(lambda (x,y):\
                cmic[x,y],list(itertools.product(mac[j],mac[i]))))
    return cmac

def update_count(cmic,cmac,mac,im,jm):
    # update macro-state count matrix
    n = len(cmic)
    m = len(mac)
    cmac_new = cmac
    for i in range(m):
        for j in [im,jm]:
            # update number of transitions FROM macrostates involved
            cmac_new[j,i] = reduce(lambda x,y: x+y,map(lambda(x,y):\
                cmic[x,y],list(itertools.product(mac[j],mac[i]))))
            # update number of transitions TO macrostates involved
            cmac_new[i,j] = reduce(lambda x,y: x+y,map(lambda(x,y):\
                cmic[x,y],list(itertools.product(mac[i],mac[j]))))

    return cmac_new

def sort_vecs(v,lw,rw):
    indx_eig= np.argsort(-np.real(v))
    v_sort = v[indx_eig]
    lw_sort = lw[:,indx_eig]
    rw_sort = rw[:,indx_eig]
    return v_sort,lw_sort,rw_sort

def test_sign(v):
    """check whether positive and negative signs are present in vector"""
    test = False
    if any(v > 0.) and  any(v<0):
        test = True 
    return test

def split_sign(macro, lvec):
    """ split based on sign structure """
    # calculate spread in eigenvector
    nt = len(macro)
    spread = []
    vals = lvec
    for k, v in macro.iteritems():
        # check that there are positive and negative values in evec
        if test_sign(vals[v]):
            #spread.append(np.sum(vals**2))
            spread.append(np.mean(vals[v]**2))
        else:
            spread.append(0.)
    isplit = np.argsort(-np.array(spread))[0]
#    print "         macrostate to split: %i"%isplit,np.array(spread)
    # split
    lvec_split = lvec[macro[isplit]]
#    print lvec_split
    elems = []
    for i in filter(lambda x: lvec_split[x] < 0.,\
        range(len(macro[isplit]))):
        elems.append(macro[isplit][i])
    macro_new = copy.deepcopy(macro)
    macro_new[nt] = elems
    # update old macrostate
    for i in elems: 
        macro_new[isplit].remove(i)
    return macro_new,vals

def split_sigma(macro,lvec):
    """ split based on distribution """
    nt = len(macro)
    spread = []
    for i in macro.keys():
        spread.append(np.std(lvec[macro[i]]))
    # split macrostates with maximum spread
    isplit = np.argsort(-np.array(spread))[0]
    #print "         macrostate to split: %i"%isplit,spread[isplit]
    # split based on distribution
    elems = []
    keep = []
    val_max =  np.max(lvec[macro[isplit]])
    val_min =  np.min(lvec[macro[isplit]])
    vals = (lvec[macro[isplit]] - val_min)/(val_max - val_min)  
    for i in filter(lambda x: vals[x] < 0.5,range(len(macro[isplit]))):
        elems.append(macro[isplit][i])
    for i in filter(lambda x: vals[x] >= 0.5,range(len(macro[isplit]))):
        keep.append(macro[isplit][i])
    macro_new = copy.deepcopy(macro)
    macro_new[nt] = elems
    #print macro_new
    # update old macrostate
    for i in elems: 
        macro_new[isplit].remove(i)
    macro = copy.deepcopy(macro_new)
    return macro,vals

def metastability(T):
    return np.sum(np.diag(T))

def beta(imc,mcsteps):
    # inverse temperature for MCSA
    x = imc - 1
    a = 10./mcsteps
    temp = (1 + (math.exp(-a*x)-1.)/(1.-math.exp(-a*mcsteps))) # MCSA temperature
    try:
        beta = 1./temp
    except ZeroDivisionError:
        beta = 1e20
    return beta

def metropolis(delta):
    if delta < 0: 
        return True
    else:
        accept = False
        p = min(1.0,np.exp(-delta))
        rand = random.random()
        if (rand < p): 
            accept = True
        return accept

#def do_mc(macro, count=None, fout="mc.dat", lagt=None):
#    print "\n Optimizing the lumped MSM\n"
#    out = open(fout, "w")
#    out.write("#    iter       q \n")
#    nstates = len(count)
#    nt = len(macro) 
#    mcsteps = len(count)*1000*nt # mc steps per block
#    mcsteps_max = nt*20000 # maximum number of mc steps 
#    count_mac = map_micro2macro(count,macro)
#    K,T = calc_rate(count_mac,lagt)
#    q =  metastability(T)
#    print " initial:", q
#    print macro
#    q_opt = q
#    macro_opt = copy.deepcopy(macro)
#    cont = True
#    nmc = 0 # number of mc blocks
#    cont = True
#    reject = 0
#    while cont:
#        imc = 0 
#        out.write ("%6i %12.10f\n"%(imc + nmc*mcsteps,q))
#        while imc < mcsteps:
#            # try ramdom insertion of a microstate in a macrostate
#            imac = 0
#            jmac = 0
#            while imc < mcsteps:
#                imc +=1
#                #print "\nStep %i"%imc
#                while True:
#                    # choose microstate to move around
#                    imic = random.choice(range(nstates))
#                    #print "      swapping micro-state: %i"%imic
#                    # find macrostate it belongs to
#                    imac = int([x for x in range(nt) if imic in macro[x]][0])
#                    if len(macro[imac]) > 1:
#                        # choose final macrostate
#                        jmac = random.choice([x for x in range(nt) if x!=imac])
#                        break
#                # move microstate from i to j
#                macro_new = copy.deepcopy(macro)
#                macro_new[imac].remove(imic)
#                macro_new[jmac].append(imic)
#                # calculate transition count matrix for new mapping
#                count_mac_new = map_micro2macro(count,macro_new)
#                Kmacro_new,Tmacro_new = calc_rate(count_mac_new,lagt)
#                # calculate metastability
#                q_new = metastability(Tmacro_new)
#                #print "Q new: %g"%q_new
#                #print "temp: %g ="%temp
#                #print " imc= %g; beta = %g"%(imc,beta)
#                delta = beta(imc,mcsteps)*(q - q_new) # calculate increment (Q is a -Energy)
#                #print delta
#                if metropolis(delta):
#                    #print "ACCEPT"
#                    macro = copy.deepcopy(macro_new)
#                    count_mac = count_mac_new
#                    q = q_new
#                    if q > q_opt:
#                        q_opt = q
#                        macro_opt = copy.deepcopy(macro)
#                else:
#                    reject+=1
#                    #print " REJECT"
#
#                if imc%100==0:
#                    out.write ("%6i %12.10e %12.10e\n"%(imc + nmc*mcsteps,q,1./beta(imc,mcsteps)))
#            nmc +=1
#        cont = False    
#    print " final :",q_opt
#    print macro_opt
#    print " acceptance:",1.-float(reject)/mcsteps
#    return macro_opt
