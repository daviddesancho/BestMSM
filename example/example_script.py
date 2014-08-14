#!/usr/bin/env python

import sys
import bestmsm.trajectory as traj
import bestmsm.msm as msm

traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

msm_ala5 = msm.MasterMSM([traj_ala5])
