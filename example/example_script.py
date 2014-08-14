#!/usr/bin/env python

import sys
import bestmsm.trajectory as traj
import bestmsm.msm as msm

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

print " Filename : %s"%traj_ala5.filename
print "    time between frames : %s"%traj_ala5.dt

msm_ala5 = msm.MSM([traj_ala5])
for l in [1., 2., 5., 10., 20., 50.]:
    lagt = traj_ala5.dt*l
    N = msm_ala5.calc_count(lagt=lagt)
    T, keep_states, keep_keys = \
        msm_ala5.calc_trans(count=N, lagt=lagt)
    rate = msm_ala5.calc_rate(trans=T, lagt=lagt)
    tauT, peqT, rvecsT, lvecsT =  msm_ala5.calc_eigs(rate=rate)


    sys.exit()
#print "    length of trajectory : %g"%len(traj_ala5.time)
#print "    states : ",traj_ala5.keys
