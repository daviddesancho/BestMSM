import unittest
import bestmsm.trajectory as traj


class Ala5TimeSeriesTest(unittest.TestCase):
    def setUp(self):
        self.traj_ala5 = traj.TimeSeries("tests/ala5_32states_timeseries.dat")
        self.state_keys = [x.split()[0] for x in open("tests/states", "r").readlines()]

    def test_filename(self):
        assert self.traj_ala5.filename == "tests/ala5_32states_timeseries.dat"

    def test_data(self):
        assert isinstance(self.traj_ala5.time, list)
        assert isinstance(self.traj_ala5.states, list)
        assert len(self.traj_ala5.states) == len(self.traj_ala5.time)

    def test_dt(self):
        assert self.traj_ala5.dt == 1

    def test_keys(self):
        print self.traj_ala5.keys
        print self.state_keys
        for i in self.state_keys:
            print "testing:", i
            assert str(i) in self.traj_ala5.keys

if __name__ == "__main__":
    unittest.main()

