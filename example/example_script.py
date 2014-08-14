#!/usr/bin/env python

import sys
import bestmsm.trajectory as traj
import bestmsm.msm as msm

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

print " Filename : %s"%traj_ala5.filename
print "    time between frames : %s"%traj_ala5.dt

msm_ala5 = msm.MSM([traj_ala5])
for l in [1, 2, 5, 10, 20, 50]:
    lagt = traj_ala5.dt*l
#    count = msm_ala5.calc_count(lagt=lagt)
    trans = msm_ala5.calc_trans( lagt=lagt)
    print count
    sys.exit()
#sys.exit()
#print "    length of trajectory : %g"%len(traj_ala5.time)
#print "    states : ",traj_ala5.keys
