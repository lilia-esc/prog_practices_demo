"""Unit tests for selection parsers"""
import analysis_demo.gmx_analysis


class SelectionTest(analysis_demo.gmx_analysis.GmxAnalysis):
    def __init__(self, structfile):
        super().__init__()
        self.u = self.load_trajectory(structfile)

    def __call__(self, selection):
        agroup = self.u.select_atoms(self.selection_parser[selection])
        aidx = agroup.atoms.indices + 1
        refidx = []

        with open("data/{}.ndx".format(selection)) as f:
            for line in f:
                refidx.extend([int(i) for i in line.strip().split()])

        print(aidx)
        print(refidx)

        assert(len(aidx) == len(refidx))
        assert(sorted(aidx) == sorted(refidx))


def test_protein():
    prot_test = SelectionTest("data/SARS_CoV_2_protease.gro")
    prot_test("protein")


def test_backbone():
    prot_test = SelectionTest("data/SARS_CoV_2_protease.gro")
    prot_test("backbone")


def test_alpha_C():
    prot_test = SelectionTest("data/SARS_CoV_2_protease.gro")
    prot_test("alpha_C")


def test_heavy():
    prot_test = SelectionTest("data/SARS_CoV_2_protease.gro")
    prot_test("heavy")


def test_alkane_ua():
    poly_test = SelectionTest("data/c35.gro")
    poly_test("alkane_ua")
