#!/usr/bin/env python

import bestmsm.trajectory as traj
import bestmsm.msm as msm
import bestmsm.pcca as pcca
import bestmsm.visual_lib as visual_lib


traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

msm_ala5 = msm.MasterMSM([traj_ala5, traj_ala5])

msm_ala5.chapman_kolmogorov(N=5)

msm_ala5.do_msm(lagt=1)

evals = msm_ala5.msms[1].tauT

#msmpcca = pcca.PCCA(msm_ala5.msms[1], N=4)

#msmpcca.optim()
