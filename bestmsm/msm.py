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

