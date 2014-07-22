#!/usr/bin/env python

import bestmsm.trajectory as traj
import bestmsm.msm as msm

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")
print " Filename : %s"%traj_ala5.filename
print "    time between frames : %s"%traj_ala5.dt
print "    length of trajectory : %g"%len(traj_ala5.time)
print "    states : ",traj_ala5.keys

msm_ala5 = msm.MSM([traj_ala5,traj_ala5],traj_ala5.dt)
