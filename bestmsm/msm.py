""" 
================================
 MSM module from BestMSM	
================================
"""
import sys
import numpy as np
import networkx as nx

class MSM:
	def __init__(self,trajectories,lag):
		# merge state keys from all trajectories
		self.merge_trajs(trajectories)
		# calculate count matrix
		self.calc_count(trajectories,lag=1)
		# find largest strongly connected sets
		self.check_connect()
		# calculate transition matrix
		self.calc_trans()
		# calculate eigenvalues and eigenvectors
		self.calc_eigs()

	def merge_trajs(self,trajectories):
		new_keys = []
		for traj in trajectories:
			new_keys += filter(lambda x: x not in new_keys,traj.keys)
		self.keys = new_keys

	def calc_count(self,trajectories,lag):
		""" calculate transition count matrix """
		keys = self.keys
		nkeys = len(keys)
		count = np.zeros([nkeys,nkeys],int)
		for traj in trajectories:
			raw = traj.states
			nraw = len(raw)
			pstate = "NULL"
			for i in range(0,nraw-lag,lag):
				j = i+lag
				state_i = raw[i]
				state_j = raw[j]
				if state_i in keys:
					idx_i = keys.index(state_i)
				if state_j in keys:
					idx_j = keys.index(state_j)
				count[idx_j][idx_i] += 1
		self.count = count

	def check_connect(self):
		""" check connectivity of rate matrix """
		D = nx.DiGraph(self.count)
		self.keep_states = sorted(nx.strongly_connected_components(D)[0])
		self.keep_keys = map(lambda x: self.keys[x],self.keep_states)

	def calc_trans(self):
		""" calculate transition matrix """
		count = self.count
		keep_states = self.keep_states
		nkeep = len(keep_states)
		T = np.zeros([nkeep,nkeep],float)
		for i in range(nkeep):
			ni = reduce(lambda x,y:x+y, map(lambda x: count[keep_states[x]][keep_states[i]], range(nkeep)))
			for j in range(nkeep):
				T[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)
		self.T = T	

	def calc_rate(self,lagt):
		""" calculate rate matrix using a Taylor series as described in De Sancho et al
		J Chem Theor Comput (2013)"""
		nkeep = len(keep_states)
		K = self.T/lagt
		for i in range(nkeep):
			K[i][i] = -(np.sum(K[:i,i]) + np.sum(K[i+1:,i]))
		self.K = K

	def get_eigs(self):
		#evalsK,rvecsK = np.linalg.eig(K)
		self.evalsK,self.lvecsK,self.rvecsK = scipyla.eig(self.K,left=True)
		self.evalsT,self.lvecsT,self.rvecsT = scipyla.eig(self.T,left=True)
		# from K
		elistK = []
	    for i in range(nkeep):
		    elistK.append([i,np.real(evalsK[i])])
		elistK.sort(msm_lib.esort)
		ieqK,eK = elistK[0]
		peqK_sum = reduce(lambda x,y:x+y, map(lambda x: rvecsK[x,ieqK], 
			range(nkeep)))
		peqK = rvecsK[:,ieqK]/peqK_sum
		# from T
		elistT = []
		for i in range(nkeep):
			elistT.append([i,np.real(evalsT[i])])
		elistT.sort(msm_lib.esort)
		ieqT,eT = elistT[0]
		peqT_sum = reduce(lambda x,y:x+y, map(lambda x: rvecsT[x,ieqT],
			range(nkeep)))
		peqT = rvecsT[:,ieqT]/peqT_sum
