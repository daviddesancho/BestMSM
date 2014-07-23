""" 
================
 Trajectory 
================
"""

class TimeSeries:
    """Time series reader

    A class to store all information about trajectories 
    to be used for the construction of the MSM

    Parameters
    ----------
    filename : str 
        The name of the file to read. It is also used for 
        indexing.

    Attributes
    ----------
    filename : str 
        The name of the file to read. It is also used for 
        indexing.

    time : list of floats
        Time stamps for time series

    states : list
        Time series of states, whatever they are.
        
    keys : list
        The names of the states.

    """
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
        """ Finds spacing between frames """
        self.dt = self.time[1] - self.time[0]

    def find_keys(self):
        """ Finds state names """
        keys = []
        for s in self.states:
            if s not in keys:
                keys.append(s)
        # sort if all integers
        if all(isinstance(x,int) == True for x in keys):
            self.keys = keys.sort()
        else:
	        self.keys = keys
