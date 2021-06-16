import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
import MDAnalysis.analysis.align
import pandas as pd


def load_trajectory(structfile, trajfile):
    return mda.Universe(structfile, trajfile)


def rgyr(coords, masses):
    com = np.average(coords, weights=masses, axis=0)
    sq_distances = np.sum((coords - com)**2, axis=1)
    Rg = np.sqrt(np.average(sq_distances, weights=masses))
    return Rg


def radius_of_gyration(u=None, selection=None):
    times = []
    Rgs = []
    sel = u.select_atoms(selection)

    for ts in u.trajectory:
        Rgval = rgyr(sel.positions, sel.masses)
        times.append(ts.time)
        Rgs.append(Rgval)

    times = np.array(times)
    Rgs = np.array(Rgs)

    return pd.DataFrame({"t": times, "Rg": Rgs})


def plot_data(df):
    fig, ax = plt.subplots(dpi=150)
    cols = df.columns.astype(str)
    ax.plot(df[cols[0]], df[cols[1]])
    ax.set_xlabel(cols[0])
    ax.set_ylabel(cols[1])

    return fig, ax


if __name__ == "__main__":
    # u = load_trajectory("../data/c35.gro", "../data/c35.xtc")
    u = load_trajectory("../data/SARS_CoV_2_protease.gro", "../data/SARS_CoV_2_protease.xtc")

    # selection = "name C*"
    selection = "name CA"
    # selection = "name CA or name C or name N"
    # selection = "protein and not name H*"
    # selection = "protein"

    df_rg = radius_of_gyration(u, selection)

    print("Writing CSV")
    # df_rg.to_csv("out/c35.csv")
    df_rg.to_csv("out/protease_alpha_C.csv")
    # df_rg.to_csv("out/protease_backbone.csv")
    # df_rg.to_csv("out/protease_heavy.csv")
    # df_rg.to_csv("out/protease_protein.csv")

    fig, ax = plot_data(df_rg)
    print("Saving plotted data")
    # fig.savefig("out/c35.png")
    fig.savefig("out/protease_alpha_C.png")
    # fig.savefig("out/protease_backbone.png")
    # fig.savefig("out/protease_heavy.png")
    # fig.savefig("out/protease_protein.png")
