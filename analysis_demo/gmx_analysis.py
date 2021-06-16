"""Class to analyse a GROMACS trajectory"""

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.align


class GmxAnalysis:
    def __init__(self):
        self.load_trajectory()

    def load_trajectory(self):
        pass

    def _rgyr(self, coords, masses):
        pass

    def radius_of_gyration(self, selection):
        pass
