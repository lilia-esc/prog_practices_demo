import analysis_demo.gmx_analysis
import argparse


class Rg(analysis_demo.gmx_analysis.GmxAnalysis):
    def __init__(self, args):
        super().__init__()
        self.u = self.load_trajectory(args.structfile, args.trajfile)
        self.selection = self.selection_parser[args.selection]

        self.outcsv = args.ocsv
        if self.outcsv is None:
            self.outcsv = 'rg.csv'

        self.outimage = args.oimage
        if self.outimage is None:
            self.outimage = 'plot.png'

    def __call__(self):
        df_rg = self.radius_of_gyration(u=self.u, selection=self.selection)
        print("Writing CSV")
        df_rg.to_csv(self.outcsv)

        fig, ax = self.plot_data(df_rg)
        print("Saving plotted data")
        fig.savefig(self.outimage)


if __name__ == "__main__":
    # Define argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("structfile", help="GROMACS structure file.")
    parser.add_argument("trajfile", help="GROMACS trajectory file.")
    parser.add_argument("selection", help="Selection group (see code for options).")
    parser.add_argument("--ocsv", help="CSV file to save dataframe to (default=rg.csv)")
    parser.add_argument("--oimage", help="Image file to save plot to (default=plot.png)")

    # Read command line arguments
    args = parser.parse_args()

    # Run analysis
    rg_calc = Rg(args)
    rg_calc()
