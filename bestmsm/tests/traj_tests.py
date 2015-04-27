import unittest
import bestmsm.trajectory as traj


class Ala5TimeSeriesTest(unittest.TestCase):
    def setUp(self):
        traj_ala5 = traj.TimeSeries("ala5_32states_timeseries.dat")
        state_keys = open("states", "r").readlines()

    def test_filename(self):
        assert self.filename == "ala5_32states_timeseries.dat"

    def test_data(self):
        assert isinstance(self.time, list)
        assert isinstance(self.states, list)
        assert len(self.states) == len(self.time)

    def test_dt(self):
        assert self.dt == 1

    def test_keys(self):
        for i in state_keys:
            assert i in self.keys

if __name__ == "__main__":
    unittest.main()

