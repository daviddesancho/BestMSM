#!/usr/bin/env python

import bestmsm.trajectory as besttraj

trajala5 = besttraj.TimeSeries("ala5_32states_timeseries.dat")
print " Filename : %s"%trajala5.filename
print "    time between frames : %s"%trajala5.dt
print "    length of trajectory : %g"%len(trajala5.time)
print "    states : ",trajala5.keys
