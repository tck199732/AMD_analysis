# I/O
import warnings

# analysis tool
import pandas as pd
import numpy as np

# miscellaneous
import random
from datetime import datetime
random.seed(datetime.now().timestamp())

# local
from pyamd.e15190 import e15190
from pyamd.utilities import root6, histo
hist_reader = root6.HistogramReader()
dfhist = histo.histogram_handler

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

    def __init__(self, fname, hname, particle, reaction):
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
        """
        self.fname = fname
        self.hname = hname
        self.particle = e15190.Particle(particle)
        self.reaction = reaction

        with root6.TFile(str(fname), mode='READ') as file:
            try:
                hist = file[hname]
            except:
                objname = file.keys()[0]
                warnings.warn(
                    f'hist name {hname} not found in root file. We will use the first object in the file, i.e. `{objname}`')
                hist = file[objname]

            self.nentries = hist.GetEntries()
            self.integral = hist.Integral()
            self.bins = (hist.GetNbinsX(), hist.GetNbinsY())
            self.binwidths = (hist.GetXaxis().GetBinWidth(1), hist.GetYaxis().GetBinWidth(1))
            self.xrange = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
            self.yrange = (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())
            self.df = hist_reader.hist2d_to_df(hist)
    
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

        df = self.correct_coverage(xrange=xrange, yrange=yrange, **kwargs) if correct_coverage else self.df
        return dfhist.projection(df, range=yrange, cut=xrange, bins=bins, normalize=True, axis='y')

    def get_normed_rapidity(self, xrange=(0,1.5), yrange=(0,800), bins=150):
        """ Project to X-axis to get distribution of `rapidity_lab` / `beam_rapidity`, i.e. normalized `dN/dy_0`
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
        return dfhist.projection(self.df, range=xrange, cut=yrange, bins=bins, normalize=True, axis='x')

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
        df = self._apply_coverage_correction(xrange=xrange, yrange=yrange, **kwargs) if correct_coverage else self.df
        
        #alias
        mass = self.particle.mass.value

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
    
    def _apply_coverage_correction(self, xrange=(0.4, 0.6), yrange=(0., 600.), seed=0):
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

        df = dfhist.query(self.df, xrange=xrange, yrange=yrange)
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
