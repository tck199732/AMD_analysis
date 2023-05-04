import iminuit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pyamd.utilities import root6, histo

hist_reader = root6.HistogramReader()
dfhist = histo.histogram_handler

class EmissionTime:
    def __init__(self, fname, hname, **kwargs):
        """ A class to analyze histogram of emission time (y-axis), the x-axis could be P, Pt, ..., etc.
        Parameters
        ----------
        fname : str
            path to the root file
        hname : str 
            name of the histogram
        **kwargs : dict
            additional attributes about the simulation, e.g. reaction, skyrme, impact_parameter
        """
        self.hname = hname
        for (key, value) in kwargs.items():
            setattr(self, key, value)

        with root6.TFile(fname, 'READ') as f:
            hist = f[hname]
            hist.SetDirectory(0)

        self.TH2D = hist
        self.nentries = hist.GetEntries()
        self.integral = hist.Integral()
        self.bins = (hist.GetNbinsX(), hist.GetNbinsY())
        self.binwidths = (hist.GetXaxis().GetBinWidth(1), hist.GetYaxis().GetBinWidth(1))
        self.xrange = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
        self.yrange = (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())
        
        self.df = hist_reader.hist2d_to_df(hist)

    def get_average_emission_time(self, range=(0,800), bins=40, option=None):

        df = dfhist.rebin2D(
            self.df, 
            range = [range, self.yrange],
            bins = [bins, self.bins[1]],
            normalize = True,
        )
        bw = (range[1] - range[0]) / bins / self.binwidths[0]

        return dfhist.profile(
            df, 
            axis = 'y', 
            option = option, 
            scale = self.nentries / (self.integral / bw),
        )

    def get_x_spectra(self):
        pass