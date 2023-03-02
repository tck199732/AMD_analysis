# I/O
import os
import pathlib
import warnings

# analysis tool
import pandas as pd
import numpy as np

# plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# miscellaneous
import random
from copy import copy
from typing import Literal
from datetime import datetime
random.seed(datetime.now().timestamp())

# local
from pyamd import PROJECT_DIR
from pyamd.e15190 import e15190
from pyamd.utilities import root6, dataframe
df_helper = dataframe.DataFrameHelper()
hist_reader = root6.HistogramReader()

class PtRapidityLAB:
    # not used, just a record
    lab_theta_range = (30.0, 75.0)  # degree
    lab_kinergy_per_A_ranges = {  # MeV/A
        'p': (20.0, 198.0),
        'd': (15.0, 263. / 2),
        't': (12.0, 312. / 3),
        '3He': (20.0, 200.0),
        '4He': (18.0, 200.0),
        '6He': (13.0, 200.0),
        '6Li': (22.0, 200.0),
        '7Li': (22.0, 200.0),
        '8Li': (22.0, 200.0),
        '7Be': (22.0, 200.0),
        '9Be': (22.0, 200.0),
        '10Be': (22.0, 200.0),
    }

    DIR = pathlib.Path(PROJECT_DIR, '/result/spectra')
    
    def __init__(self, fname, hname, particle, reaction, **kwargs):
        """ A class for handling TH2D tranverse momentum / A v.s. rapidity / beam-rapidity
        Parameters
        ----------
        fname : str or pathlib.Path
            path to the root file
        hname : str
            name of the histogram in the ROOT file
        particle : str
            particle name, e.g. 'p', 'd', 't', '3He', '4He'
        reaction : str
            reaction name, e.g. 'Ca48Ni64E140'
        skyrme : str, optional
            skyrme name, e.g. 'SkM', 'SLy4', 'SLy4_L108', by default None
        mode : str, optional
            mode name, e.g. 'primary', 'secondary', by default None
        impact_parameter : float, optional
            impact parameter, by default None
        multiplicity : int, optional
            multiplicity, by default None
        """

        self.fname = fname
        self.hname = hname
        self.particle = e15190.Particle(particle)
        self.reaction = reaction

        attributes = ['skyrme', 'mode', 'impact_parameeter', 'multiplicity']
        for attr in attributes:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])

        with root6.TFile(str(fname), mode='READ') as file:
            try:
                hist = file[hname]
            except:
                objname = file.keys()[0]
                warnings.warn(
                    f'hist name {hname} not found in root file. We will use the first object in the file, i.e. `{objname}`')
                hist = file[objname]

            self.nentries = hist.GetEntries()
            self.bins = (hist.GetNbinsX(), hist.GetNbinsY())
            self.binwidths = (hist.GetXaxis().GetBinWidth(1), hist.GetYaxis().GetBinWidth(1))
            self.xrange = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
            self.yrange = (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())
            self.df = hist_reader.hist2d_to_df(hist)

    def _slice(self, xrange=(0.4,0.6), yrange=(0., 600.)):
        """ Slice the dataframe by the given range, 
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0.4,0.6)
        yrange : tuple, optional
            y range, by default (0., 600.)
        Returns
        -------
        pandas.DataFrame
            sliced dataframe
        Notes
        -----
        if bin width is `0.005` and we want to slice from `0.4` to `0.6`,  then 0.396 is included in ROOT Projection.
        0.395 is not included becasue any value at the interval would be filled in the bin on the left. so (> and <=)
        """
        df = self.df.copy()
        x1 = xrange[0] - self.binwidths[0] / 2.
        x2 = xrange[1] + self.binwidths[0] / 2
        y1 = yrange[0] - self.binwidths[1] / 2.
        y2 = yrange[1] + self.binwidths[1] / 2
        return df.query(f'x > {x1} & x <= {x2} & y > {y1} & y <= {y2}')

    def correct_coverage(self, xrange=(0.4, 0.6), yrange=(0., 600.), seed=0):
        """ Apply coverage-correction to the histogram. 
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0.4, 0.6)
        yrange : tuple, optional
            y range, by default (0., 600.)
        seed : int, optional
            random seed, by default 0
        Returns
        ------- 
        pandas.DataFrame
            corrected dataframe
        """
        rng = np.random.default_rng(seed=seed)

        df = self._slice(xrange, yrange)
        for _, subdf in df.groupby('y'):
            if subdf.shape[0] == 0:
                df['z'].iloc[subdf.index] = 0.
                df['z_err'].iloc[subdf.index] = 0.
            else:   
                x1 = np.min(subdf['x']) + self.binwidths[0] * rng.uniform(-0.5, 0.5)
                x2 = np.max(subdf['x']) + self.binwidths[0] * rng.uniform(-0.5, 0.5)
                df['z'] /= abs(x2-x1) / np.diff(xrange)
                df['z_err'] /= abs(x2-x1) / np.diff(xrange)

        return df

    def PtSpectrum(self, xrange=(0.4, 0.6), yrange=(0., 600.), bins=30,  correct_coverage=False, **kwargs):
        return self.get_pt(xrange, yrange, bins, correct_coverage, **kwargs)

    def get_pt(self, xrange=(0.4, 0.6), yrange=(0., 600.), bins=30,  correct_coverage=False, **kwargs):
        """ Get the transverse momentum spectrum
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0.4, 0.6)
        yrange : tuple, optional
            y range, by default (0., 600.)
        bins : int, optional
            number of bins, by default 30
        correct_coverage : bool, optional   
            apply coverage correction, by default False
        **kwargs : dict
            keyword arguments for `correct_coverage`
        Returns
        -------
        pandas.DataFrame
            transverse momentum spectrum
        """

        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=yrange, **kwargs)
        else:
            df = self._slice(xrange, yrange)
        
        hist, edges = np.histogram(df['y'], range=yrange, bins=bins, weights=df['z'])
        histerr, _ = np.histogram(df['y'], range=yrange, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)

        norm = np.diff(yrange) / bins * np.diff(xrange)

        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist / norm,
            'y_err': histerr / norm,
            'y_ferr':  np.divide(histerr, hist, where=(hist != 0.0))
        })

    def RapidityDistribution(self, xrange=(0,1), yrange=(0,800), bins=100):
        return self.get_normed_rapidity(xrange, yrange, bins)
    
    def get_normed_rapidity(self, xrange=(0,1), yrange=(0,800), bins=100):
        """ Project to X-axis to get distribution of `rapidity_lab` / `beam_rapidity`
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0,1)
        yrange : tuple, optional
            y range, by default (0,800)
        bins : int, optional
            number of bins, by default 100
        Returns
        -------
        pandas.DataFrame
            rapidity distribution
        """
        df = self._slice(xrange, yrange)
        hist, _ = np.histogram(df['x'], weights=df['z'], bins=bins, range=xrange)
        hist_err, _ = np.histogram(df['x'], weights=df['z_err']**2, bins=bins, range=xrange)
        hist_err = np.sqrt(hist_err)
        
        x = np.linspace(*xrange, bins+1)
        return pd.DataFrame({
            'x' : 0.5 * (x[1:] + x[:-1]),
            'y' : hist,
            'y_err' : hist_err,
            'y_ferr' : np.divide(hist_err, hist, where=(hist != 0.))
        })

    def get_kinergy_theta_CMS(self, xrange=(0.4, 0.6), yrange=(0., 600.), correct_coverage=False, **kwargs):
        """ Get the 2D histogram kinergy per A versus theta in CMS
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0.4, 0.6)
        yrange : tuple, optional
            y range, by default (0., 600.)
        correct_coverage : bool, optional
            apply coverage correction, by default False
        **kwargs : dict
            keyword arguments for `correct_coverage`
        Returns
        -------
        pandas.DataFrame
            kinergy per A versus theta in CMS
        """
        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=yrange, **kwargs)
        else:
            df = self._slice(xrange, yrange)
        
        #alias
        mass = self.particle.mass

        # randomize each bin
        rng = np.random.default_rng(seed=kwargs.get('seed', 0))
        x = df['x'].to_numpy() + self.binwidths[0] * rng.uniform(-0.5, 0.5, size=df.shape[0])
        y = df['y'].to_numpy() + self.binwidths[1] * rng.uniform(-0.5, 0.5, size=df.shape[0])

        # get reaction constants
        beam_rapidity = e15190.CollisionReaction.get_rapidity_beam(self.reaction)
        beta_cms = e15190.CollisionReaction.get_betacms(self.reaction)

        # calculate rapidity in CMS and total transverse momentum
        rapidity_cms = x * beam_rapidity
        rapidity_cms -= 0.5 * np.log((1 + beta_cms) / (1 - beta_cms))
        pt = y * (self.particle.N + self.particle.Z)
        
        # calculate pz and theta in CMS
        pz_cms = pt * np.sinh(rapidity_cms)
        theta_cms = np.degrees(np.arctan2(pt, pz_cms))

        # calculate kinergy (per A) in CMS
        kinergy_cms = np.sqrt(pt**2 + mass**2) * np.cosh(rapidity_cms) - mass
        kinergy_cms_per_A = kinergy_cms / (self.particle.N + self.particle.Z)
        
        return pd.DataFrame({
            'x': theta_cms,
            'y': kinergy_cms_per_A,
            'z': df['z'],
            'z_err': df['z_err'],
            'z_ferr': df['z_ferr']
        })

    def get_pmag_rapidity_CMS(self, xrange=(0.4, 0.6), yrange=(0., 600.), correct_coverage=False, **kwargs):
        """ Get the 2D histogram momentum magnitude versus rapidity in CMS
        Parameters
        ----------
        xrange : tuple, optional
            x range, by default (0.4, 0.6)
        yrange : tuple, optional
            y range, by default (0., 600.)
        correct_coverage : bool, optional       
            apply coverage correction, by default False
        **kwargs : dict 
            keyword arguments for `correct_coverage`
        Returns
        -------
        pandas.DataFrame    
            momentum magnitude versus rapidity in CMS
        """
        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=yrange, **kwargs)
        else:
            df = self._slice(xrange, yrange)
        
        #alias
        mass = self.particle.mass

        # randomize each bin
        rng = np.random.default_rng(seed=kwargs.get('seed', 0))
        x = df['x'].to_numpy() + self.binwidths[0] * rng.uniform(-0.5, 0.5, size=df.shape[0])
        y = df['y'].to_numpy() + self.binwidths[1] * rng.uniform(-0.5, 0.5, size=df.shape[0])

        # get reaction constants
        beam_rapidity = e15190.CollisionReaction.get_rapidity_beam(self.reaction)
        beta_cms = e15190.CollisionReaction.get_betacms(self.reaction)

        rapidity_cms = x * beam_rapidity
        rapidity_cms -= 0.5 * np.log((1 + beta_cms) / (1 - beta_cms))
        pt = y * (self.particle.N + self.particle.Z)
        
        # calculate pz and theta in CMS
        pz_cms = pt * np.sinh(rapidity_cms)
        pmag_cms = np.sqrt(pt**2 + pz_cms**2)
        pmag_cms_per_A = pmag_cms / (self.particle.N + self.particle.Z)

        return pd.DataFrame({
            'x': rapidity_cms,
            'y': pmag_cms_per_A,
            'z': df['z'],
            'z_err': df['z_err'],
            'z_ferr': df['z_ferr']
        })
        



