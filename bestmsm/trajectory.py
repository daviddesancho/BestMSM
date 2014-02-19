import sys,os

""" 
================================
 Trajectory module from BestMSM	
================================
 Functionalities:
"""

class TimeSeries:
	def __init__(self,filename):
		self.filename = filename
		self.readtraj()

	def readtraj(self):
	"""	Reads trajectory files assuming two column file corresponding to time
	and state. If there are not two columns then states are read from the
	first column and times are line numbers"""
	
