import unittest
import numpy as np
import bestmsm.msm as msm
import bestmsm.pcca as pcca


class TestTwoClusters(unittest.TestCase):
    def setUp(self):
        """
        We set up a test involving a fake transition count matrix

        """

        data = None
        self.keys = range(6)
        lagt = 1
        self.msmS = msm.MSM(data, keys=self.keys, lagt=1)
        self.msmS.count = np.array([[ 6000, 300, 0, 0, 0, 0], \
                     [ 300, 1000, 100, 0, 0, 0], \
                     [ 0, 100, 1000, 3, 0, 0], \
                     [ 0, 0, 3, 1000, 100, 0], \
                     [ 0, 0, 0, 100, 1000, 300], \
                     [ 0, 0, 0, 0, 300, 90000]])
        self.msmS.keep_states, self.msmS.keep_keys = self.msmS.check_connect()
        self.msmS.do_trans()

    def test_eigen(self):
        msmpcca = pcca.PCCA(self.msmS, N=2)
        print msmpcca.macros
        i = [x for x in msmpcca.macros.keys() if 0 in msmpcca.macros[x]][0]
        assert 1 in  msmpcca.macros[i]
        assert 2 in  msmpcca.macros[i]
        j = [x for x in msmpcca.macros.keys() if 3 in msmpcca.macros[x]][0]
        assert 4 in  msmpcca.macros[j]
        assert 5 in  msmpcca.macros[j]

if __name__ == "__main__":
    unittest.main()
