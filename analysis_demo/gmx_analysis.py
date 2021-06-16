"""Class to analyse a GROMACS trajectory"""

import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
import MDAnalysis.analysis.align
import pandas as pd


class GmxAnalysis:
    def __init__(self):
        # Class variables to be defined by child class:
        self.u = None
        self.selection_parser = {
            'protein': """protein""",
            'backbone': """name CA or name C or name N""",
            'alpha_C': """name CA""",
            'heavy': "protein and not name H*",
            'alkane_ua': """name C*"""
        }

    def load_trajectory(self, structfile, trajfile):
        """
        Loads trajectory from the given structure and trajectory files.

        Args:
            structfile (str): Path to GROMACS structure file (.gro or .tpr).
            structfile (str): Path to GROMACS trajectory file (.trr or .xtc).
        """
        return mda.Universe(structfile, trajfile)

    def _rgyr(self, coords, masses):
        r"""
        Calculates the radius of gyration of atoms with given coordinates and
        masses.

        If $r_i(t)$ is the selection of coordinates,
        $$R_g = \sqrt{\frac{1}{N} \sum_{k=1}^{N} (r_k - r_{com})^2 }$$

        Args:
            coords (np.array): Array of shape (N, 3) containing atomic coordinates.
            masses (np.array): Array of shape (N,) containing atomic masses.

        Returns:
            Radius of gyration (np.float).

        Raises:
            ValueError if coords does not have the right shape, or if the
            lengths of coords and masses do not match.
        """
        com = np.average(coords, weights=masses, axis=0)
        sq_distances = np.sum((coords - com)**2, axis=1)
        Rg = np.sqrt(np.average(sq_distances, weights=masses))
        return Rg

    def radius_of_gyration(self, u=None, selection=None):
        """
        Calculates radius of gyration of selection along trajectory.

        Args:
            u (mda.Universe): Universe (default=self.u, i.e. main universe).
            selection (str): MDAnalysis selection string (default=None, i.e. all atoms).

        Returns:
            Time-indexed radii of gyration (pd.DataFrame)
        """
        if u is None:
            u = self.u

        times = []
        Rgs = []
        sel = u.select_atoms(selection)

        for ts in u.trajectory:
            Rgval = self._rgyr(sel.positions, sel.masses)
            times.append(ts.time)
            Rgs.append(Rgval)

        times = np.array(times)
        Rgs = np.array(Rgs)

        return pd.DataFrame({"Time (ps)": times, "Rg (A)": Rgs})

    def plot_data(self, df):
        """Plots data from dataframe

        Args:
            df: Dataframe to plot. Index will be plotted on the

        Returns:
            tuple(fig, ax):
                fig: Matplotlib figure object
                ax: Matplotlib axis object corresponding to the figure
        """
        fig, ax = plt.subplots(dpi=150)
        cols = df.columns.astype(str)
        ax.plot(df[cols[0]], df[cols[1]])
        ax.set_xlabel(cols[0])
        ax.set_ylabel(cols[1])

        return fig, ax
