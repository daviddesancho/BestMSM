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

    dt : float
        Lapse between snapshots

    """
    def __init__(self, filename=None, keys=None):
        if filename:
            self.filename = filename
            self.time, self.states = self.read_traj()
            self.dt = self.find_dt()
            if keys is None:
                self.keys = self.find_keys()
            else:
                self.keys = filter(lambda x: x[0] not in ['#', '@'],
                         open(keys).readlines())
        else: 
            # asume overriding initialization
            pass

    def read_traj(self):
        """	Reads trajectory files assuming two columns: time and state. 
        If there is only a single column, times are assumed to be the line
        numbers

        """
        raw = filter(lambda x: x[0] not in ['#', '@'],
                     open(self.filename).readlines())
        try:
            data = map(lambda x: (x.split()[0], x.split()[1]), raw)
            time = [float(x[0]) for x in data]
            states = [x[1] for x in data]
        except IndexError:
            data = map(lambda x: x.split()[0], raw)
            states = [x for x in data]
            time = range(len(states))
        return time, states
 
    def find_dt(self):
        """ Finds spacing between frames """
        return self.time[1] - self.time[0]

    def find_keys(self):
        """ Finds state names """
        keys = []
        for s in self.states:
            if s not in keys:
                keys.append(s)
        # sort if all integers
        if all(isinstance(x,int) is True for x in keys):
            keys = keys.sort()
        else:
            pass
        return keys
