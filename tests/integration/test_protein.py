"""Integration test for protein calculation"""
from analysis_demo import gmx_analysis
import numpy as np


class ProteinTest(gmx_analysis.GmxAnalysis):
    def __init__(self):
        super().__init__()
        self.u = self.load_trajectory("data/SARS_CoV_2_protease.gro", "data/SARS_CoV_2_protease.xtc")
        self.selection_alpha_C = self.selection_parser["alpha_C"]
        self.selection_backbone = self.selection_parser["backbone"]
        self.selection_heavy = self.selection_parser["heavy"]
        self.selection_protein = self.selection_parser["protein"]

    def compare_alpha_C(self):
        """Compares outputs of GmxAnalysis Rg of protein C-alpha atoms
        to g_gyrate outputs within tolerance of 1e-4."""
        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection_alpha_C)

        times = df_rg.t.to_numpy()
        print(times)  # print logging
        rg = df_rg.Rg.to_numpy()
        print(rg)   # print logging

        gmx_times = []
        gmx_rg = []

        f = open("data/protease_Rg_alpha-C.xvg")
        for line in f:
            if line.strip()[0] not in ['#', '@']:
                gmx_times.append(float(line.strip().split()[0]))
                gmx_rg.append(10 * float(line.strip().split()[1]))
        f.close()

        print(times - gmx_times)  # print logging
        print(rg - gmx_rg)  # print logging

        assert(np.allclose(times, gmx_times, atol=1e-4))
        assert(np.allclose(rg, gmx_rg, atol=1e-4))

    def compare_backbone(self):
        """Compares outputs of GmxAnalysis Rg of protein backbone atoms
        to g_gyrate outputs within tolerance of 1e-2."""
        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection_backbone)

        times = df_rg.t.to_numpy()
        print(times)  # print logging
        rg = df_rg.Rg.to_numpy()
        print(rg)   # print logging

        gmx_times = []
        gmx_rg = []

        f = open("data/protease_Rg_backbone.xvg")
        for line in f:
            if line.strip()[0] not in ['#', '@']:
                gmx_times.append(float(line.strip().split()[0]))
                gmx_rg.append(10 * float(line.strip().split()[1]))
        f.close()

        print(times - gmx_times)  # print logging
        print(rg - gmx_rg)  # print logging

        assert(np.allclose(times, gmx_times, atol=1e-2))  # differences in mass
        assert(np.allclose(rg, gmx_rg, atol=1e-2))  # differences in mass

    def compare_heavy(self):
        """Compares outputs of GmxAnalysis Rg of heavy atoms
        to g_gyrate outputs within tolerance of 1e-1."""
        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection_heavy)

        times = df_rg.t.to_numpy()
        print(times)  # print logging
        rg = df_rg.Rg.to_numpy()
        print(rg)   # print logging

        gmx_times = []
        gmx_rg = []

        f = open("data/protease_Rg_heavy_GMX.xvg")
        for line in f:
            if line.strip()[0] not in ['#', '@']:
                gmx_times.append(float(line.strip().split()[0]))
                gmx_rg.append(10 * float(line.strip().split()[1]))
        f.close()

        print(times - gmx_times)  # print logging
        print(rg - gmx_rg)  # print logging

        assert(np.allclose(times, gmx_times, atol=1e-1))  # differences in mass
        assert(np.allclose(rg, gmx_rg, atol=1e-1))  # differences in mass

    def compare_protein(self):
        """Compares outputs of GmxAnalysis Rg of protein
        to g_gyrate outputs within tolerance of 2e-1."""
        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection_protein)

        times = df_rg.t.to_numpy()
        print(times)  # print logging
        rg = df_rg.Rg.to_numpy()
        print(rg)   # print logging

        gmx_times = []
        gmx_rg = []

        f = open("data/protease_Rg_all_GMX.xvg")
        for line in f:
            if line.strip()[0] not in ['#', '@']:
                gmx_times.append(float(line.strip().split()[0]))
                gmx_rg.append(10 * float(line.strip().split()[1]))
        f.close()

        print(times - gmx_times)  # print logging
        print(rg - gmx_rg)  # print logging

        assert(np.allclose(times, gmx_times, atol=2e-1))  # differences in mass
        assert(np.allclose(rg, gmx_rg, atol=2e-1))  # differences in mass


def test_alpha_C():
    test = ProteinTest()
    test.compare_alpha_C()


def test_backbone():
    test = ProteinTest()
    test.compare_backbone()


def test_heavy():
    test = ProteinTest()
    test.compare_heavy()


def test_protein():
    test = ProteinTest()
    test.compare_protein()
