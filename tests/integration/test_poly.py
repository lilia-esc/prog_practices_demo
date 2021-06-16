"""Integration test for polymer calculation"""
from analysis_demo import gmx_analysis
import numpy as np


class PolyTest(gmx_analysis.GmxAnalysis):
    def __init__(self):
        super().__init__()
        self.u = self.load_trajectory("data/c35.gro", "data/c35.xtc")
        self.selection = self.selection_parser["alkane_ua"]

    def __call__(self):
        """Compares output of GmxAnalysis Rg
        to g_gyrate output within tolerance of 1e-5."""

        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection)

        times = df_rg.t.to_numpy()
        print(times)  # print logging
        rg = df_rg.Rg.to_numpy()
        print(rg)   # print logging

        gmx_times = []
        gmx_rg = []

        f = open("data/c35_Rg_GMX.xvg")
        for line in f:
            if line.strip()[0] not in ['#', '@']:
                gmx_times.append(float(line.strip().split()[0]))
                gmx_rg.append(10 * float(line.strip().split()[1]))
        f.close()

        print(times - gmx_times)  # print logging
        print(rg - gmx_rg)  # print logging

        assert(np.allclose(times, gmx_times, atol=1e-5))
        assert(np.allclose(rg, gmx_rg, atol=1e-5))


def test_poly():
    test = PolyTest()
    test()
