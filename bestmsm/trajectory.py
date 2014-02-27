""" 
================================
 Trajectory module from BestMSM	
================================
"""

class TimeSeries:
    def __init__(self, filename, keys=None):
        self.filename = filename
        self.readtraj()
        self.find_dt()
        if keys == None:
            self.find_keys()
        else:
            self.keys = filter(lambda x: x[0] not in ['#','@'],
        	    open(keys).readlines())
    
    def readtraj(self):
        """	Reads trajectory files assuming two columns: time and state. 
        If there is only a single column, times are assumed to be the line
        numbers """
        raw = filter(lambda x: x[0] not in ['#','@'],
            open(self.filename).readlines())
        try:
            data = map(lambda x: (x.split()[0],x.split()[1]), raw)
            time = [float(x[0]) for x in data]
            states = [x[1] for x in data]
        except IndexError:
            states = [x[0] for x in data]
            time = range(len(state))
        self.time = time
        self.states = states
 
    def find_dt(self):
        """ finds spacing between frames """
        self.dt = self.time[1] - self.time[0]

    def find_keys(self):
        """ finds state names """
        keys = []
        for s in self.states:
            if s not in keys:
                keys.append(s)
        # sort if all integers
        if all(isinstance(x,int) == True for x in keys):
            self.keys = keys.sort()
        else:
	        self.keys = keys
