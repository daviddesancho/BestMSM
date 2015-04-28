import unittest
import numpy as np
import bestmsm.msm as msm


class TestSinglePathway(unittest.TestCase):
    def setUp(self):
        """
        We set up a test involving a fake transition count matrix

        """

        data = None
        self.keys = range(6)
        lagt = 1
        self.msmS = msm.MSM(data, keys=self.keys, lagt=1)
        self.msmS.count = np.array([[ 6000, 3, 0, 0, 0, 0], \
                     [ 3, 1000, 3, 0, 0, 0], \
                     [ 0, 3, 1000, 3, 0, 0], \
                     [ 0, 0, 3, 1000, 3, 0], \
                     [ 0, 0, 0, 3, 1000, 3], \
                     [ 0, 0, 0, 0, 3, 90000]])

    def test_connect(self): 
        self.msmS.keep_states, self.msmS.keep_keys = self.msmS.check_connect()
        assert [x in self.msmS.keep_states for x in self.keys]
        assert [x in self.msmS.keep_keys for x in self.keys]

    def test_trans(self):
        self.msmS.keep_states, self.msmS.keep_keys = self.msmS.check_connect()
        self.msmS.do_trans()
        for pi in self.msmS.peqT:
            assert pi > 0.
#CP 
#6;000 2 2 0 0 0
#2 1;000 0 2 2 0
#2 0 1;000 2 2 0
#0 2 2 1;000 0 2
#0 2 2 0 1;000 2
#0 0 0 2 2 90;000
#2
#666666664
#3
#777777775
#CH 
#3;500 2 0 0 0 0
#2 1;000 2 0 0 2
#0 2 1;000 2 0 2
#0 0 2 1;000 2 2
#0 0 0 2 3;500 0
#0 2 2 2 0 90;000

if __name__ == "__main__":
    unittest.main()
    
