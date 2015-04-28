import unittest
import numpy as np
import bestmsm.msm as msm
import bestmsm.pcca as pcca
import bestmsm.trajectory as traj
import bestmsm.pcca_lib as pcca_lib


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

class TestDihedral(unittest.TestCase):
    def setUp(self):
        """
        Test involving a Ramachandran angle shifting from A to E

        """
        # read trajectory
        self.traj2states = traj.TimeSeries(filename="tests/alaTB_traj_toy.dat")
        data = [self.traj2states]
        # generate MSM
        self.msm2 = msm.MSM(data, keys=self.traj2states.keys, lagt=10)
        self.msm2.do_count()
        self.msm2.do_trans()

    def test_mapping(self):
        # do Perron clusters
        msm2pcca = pcca.PCCA(self.msm2, N=2)
        # find clusters based on putative cluster centers
        mac1 = [x for x in msm2pcca.macros.keys() if self.msm2.keep_keys.index('1') in msm2pcca.macros[x]][0]
        mac2 = [x for x in msm2pcca.macros.keys() if self.msm2.keep_keys.index('4') in msm2pcca.macros[x]][0]
        # check that everything is in the right place
        for i in range(1,2):
            assert self.msm2.keep_keys.index(str(i)) in msm2pcca.macros[mac1] 
        for i in range(4,5):
            assert self.msm2.keep_keys.index(str(i)) in msm2pcca.macros[mac2]
        # map trajectory onto clusters
        msm2pcca.map_trajectory()
        msm2pcca.do_count()
        # map transition count matrix onto clusters
        recount = pcca_lib.map_micro2macro(self.msm2.count, msm2pcca.macros, self.msm2.keep_states)
        assert (msm2pcca.count == recount).all()

if __name__ == "__main__":
    unittest.main()
