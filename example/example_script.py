#!/usr/bin/env python

import bestmsm.trajectory as traj
import bestmsm.msm as msm
import bestmsm.pcca as pcca

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

msm_ala5 = msm.MasterMSM([traj_ala5, traj_ala5])

msm_ala5.do_msm(lagt=1)
msm_ala5.msms[1]
msm_ala5.msms[1].do_count()
msm_ala5.msms[1].do_trans()
msm_ala5.msms[1].boots()

#time,popul = msm_ala5.msms[1].propagateK(init=['11111'])
#time,popul = msm_ala5.msms[1].propagateK(p0="pini.dat")

#msmpcca = pcca.PCCA(msm_ala5.msms[1], N=4)

#msmpcca.optim()
