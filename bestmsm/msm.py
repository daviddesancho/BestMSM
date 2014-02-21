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
		self.keep_keys = map(lambda x: self.keys[x],keep_states)

	def calc_trans(keep_states,state_keys,count):
		""" calculate transition matrix """
		nkeep = len(keep_states)
		T = np.zeros([nkeep,nkeep],float)
		n_nonzero = 0
		neighbours = {}
		for i in range(nkeep):
			ni = reduce(lambda x,y:x+y, map(lambda x: count[keep_states[x]][keep_states[i]], range(nkeep)))
			for j in range(nkeep):
				T[j][i] = float(count[keep_states[j]][keep_states[i]])/float(ni)
				#if count[keep_states[j]][keep_states[i]] != 0:
				#	n_nonzero += 1
				#	diff = difference(state_keys[keep_states[i]],state_keys[keep_states[j]])
				#	if diff not in neighbours.keys():
				#		neighbours[diff] = 0
				#	neighbours[diff] += 1
		#print "  number of nonzero elements in transition matrix = %i"%n_nonzero
		#print "                  total size of transition matrix = %i"%(nkeep**2)
	
		return T,neighbours

