#!/usr/bin/env python

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from subprocess import call
from multiprocessing import Process

beta = 1./(300*8.314e-3)

class FourState:
    """ numerical example of a four state model"""
    
    def __init__(self,option=None):
        if option == "random":
                self.ran_initialize()
        else:
                self.fixed_initialize()
        self.calc_trans()
        self.calc_rate()
        self.calc_eigs()
    
    def ran_initialize(self):
        """ initialize count matrix """
        # initialize randomly
        randint = np.random.randint(1000,size=(4,4)) 
        # disconnect end states
        randint[0,3] = 0 ; randint[3,0] = 0
        # enforce detailed balance
        randint = np.triu(randint)
        for i in range(1,4):
                for j in range(i):
                        randint[i,j] = randint[j,i]
        # enforce metastability of states
        for i in range(4): 
                randint[i,i] = randint[i,i]**2
        diagsort = np.sort(np.diag(randint)) # enforce FF and UU more stable than I1 and I2
        randint[0,0] = diagsort[2]
        randint[1,1] = diagsort[0]
        randint[2,2] = diagsort[1]
        randint[3,3] = diagsort[3]
        
        self.count = randint
        print "\n count matrix:"
        print self.count
    
    def fixed_initialize(self):
        """ initialize count matrix """
        # initialize randomly
        self.count = np.array([[468324  ,  832 ,   470 ,     0],\
                        [   832,  120096 ,   350  ,   52],\
                        [   470 ,  350, 92721 ,   1831],\
                        [     0  ,  52 ,   1831, 636804]])
        print "\n count matrix:"
        print self.count
    
    
    def calc_trans(self):
        """ calculate transition matrix (assume lag time is 1)"""
        count = self.count
        trans = np.zeros((4,4),float)
        for i in range(4):
            nsum = np.sum(count[:,i])
            for j in range(4):
                trans[j,i] = float(self.count[j,i])/nsum
        self.trans = trans
        print "\n transition matrix:"
        print self.trans
    
    def calc_rate(self):
        """ calculate rates """
        K = self.trans
        for i in range(4):
                K[i,i] = 0.
                K[i,i] = -np.sum(K[:,i])
        self.K = K
        print "\n rate:"
        print self.K
    
    def calc_eigs(self):
            """ calculate eigenvectors """
            evals,rvecs = np.linalg.eig(self.K)
            lvecs = np.linalg.inv(rvecs)
            # sort eigenvectors
            order = np.argsort(-evals)
            evals_sort = evals[order]
            rvecs_sort = np.zeros((4,4),float)
            for i in range(4):
                    rvecs_sort[:,i] = rvecs[:,order[i]]
            lvecs_sort = np.linalg.inv(rvecs_sort)
            self.peq = rvecs_sort[:,0]/np.sum(rvecs_sort[:,0])
            print "\n equilibrium probabilities:"
            print self.peq
            self.evals = evals_sort
            self.rvecs = rvecs_sort
            self.lvecs = lvecs_sort
            print "\n eigenvalues:"
            print self.evals
    
    def run_simul(self,label=None):
            """ simulate a relaxation from unfolded to folded"""
            p0 = [1.,0.,0.,0.]
            pt = []
            logt = np.arange(0,6,0.1)
            time = 10**logt
            for t in time:
                    expdiagkt = np.diag(np.exp(self.evals*t))
                    expKt = np.dot(self.rvecs,np.dot(expdiagkt,self.lvecs))
                    pt.append(np.dot(expKt,p0))
            data = []
            for i in range(4):
                    data.append([time,map(lambda x: x[i],pt)])
            p = Process(target=plot_graph, args=[data,label])
            p.start()
    
    def run_commit(self):
            """ calculate committors and reactive flux """
            K = self.K
            peq = self.peq
            # define end-states
            UU = [0]
            FF = [3]
            UUFF = UU+FF
            I = filter(lambda x: x not in UU+FF, range(4))
            NI = len(I)
    
            # calculate committors
            b = np.zeros([NI], float)
            A = np.zeros([NI,NI], float)
            for j_ind in range(NI):
                    j = I[j_ind]
                    sum = 0.
                    for i in FF:
                            sum+= K[i][j]
                    b[j_ind] = -sum
                    for i_ind in range(NI):
                            i = I[i_ind]
                            A[j_ind][i_ind] = K[i][j]               
            # solve Ax=b
            Ainv = np.linalg.inv(A)
            x = np.dot(Ainv,b)
            XX = np.dot(Ainv,A)
    
            pfold = np.zeros(4,float)
            for i in range(4):
                    if i in UU:
                            pfold[i] = 0.0
                    elif i in FF:
                            pfold[i] = 1.0
                    else:
                            ii = I.index(i)
                            pfold[i] = x[ii]
            self.pfold = pfold
            print "\n pfold values:"
            print pfold
                                                    
            # stationary distribution
            pss = np.zeros([4],float)
            for i in range(4):
                    pss[i] = (1-pfold[i])*peq[i]
            self.pss = pss
    
            # flux matrix and reactive flux
            J = np.zeros([4,4],float)
            sum_flux = 0
            for i in range(4):
                for j in range(4):
                            J[j][i] = K[j][i]*peq[i]*(pfold[j]-pfold[i])
                            if j==3: #  dividing line corresponds to I to F transitions
                                    sum_flux += J[j][i]
    
            print "\n flux :"
            print J
            print "\n reactive flux: %g"%sum_flux
            self.J = J
            self.sum_flux = sum_flux
