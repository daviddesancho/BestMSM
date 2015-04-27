import unittest
import bestmsm.trajectory as traj


class Ala5TimeSeriesTest():
    def setUp(self):
        traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")

    def test_filename(self):
        assert self.filename == "ala5_32states_timeseries.dat"

    def test_data(self):
        assert isinstance(self.dt, list)
        assert isinstance(self.states, list)

    def test_dt(self):
        assert self.dt > 0

if __name__ == "__main__":
