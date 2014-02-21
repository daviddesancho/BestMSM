""" 
================================
 MSM module from BestMSM	
================================
"""
import sys
import numpy as np

class MSM:
	def __init__(self,trajectories,lag):
		self.merge_trajs(trajectories)
		self.count = self.calc_count(trajectories,lag=1)

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
