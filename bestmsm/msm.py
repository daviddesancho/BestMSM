""" 
================================
 MSM module from BestMSM	
================================
"""
import numpy as np

class MSM:
	def __init__(self,trajectories):
		self.count = self.calc_count(trajectories,lag=1)
#		process_trajectories(trajectories)

	def calc_count(self,trajectories,lag):
		""" calculate transition count matrix """
		for traj in trajectories:
			nstates = len(traj.keys)
			count = np.zeros([nstates,nstates],int)
			nraw = len(traj.time)
			raw = traj.assign
			pstate = "NULL"
			for i in range(0,nraw-lag,lag):
				j = i+lag
				idx_i = raw[i]
				idx_j = raw[j]
#				if state_i in keys:
#					idx_i = keys.index(state_i)
#				if state_j in keys:
#					idx_j = keys.index(state_j)
#				try:
				count[idx_j][idx_i] += 1
#				except UnboundLocalError:
#					pass
#					#print " Wrong assignment"
#					#print raw[i],'-->',raw[j]
		return count
