""" 
================================
 Trajectory module from BestMSM	
================================
"""

class TimeSeries:
	def __init__(self,filename):
		self.filename = filename
		self.readtraj()

	def readtraj(self):
		"""	Reads trajectory files assuming two columns: time and state. 
		If there is only a single column, times are assumed to be the line
		numbers """
		raw = filter(lambda x: x[0] not in ['#','@'],
				open(self.filename).readlines())
		try:
			data = map(lambda x: (x.split()[0],x.split()[1]),raw)
			time = [float(x[0]) for x in data]
			states = [float(x[1]) for x in data]
		except IndexError:
			states = [float(x[0]) for x in data]
			time = range(len(state))
		self.time = time
		self.states = states
		
	def find_dt(self):
		""" finds spacing between frames """

