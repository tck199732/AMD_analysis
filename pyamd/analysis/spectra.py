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
    
    def __init__(self, particle, path=None, reaction='Ca48Ni64E140', skyrme=None, mode:Literal['primary','secondary', 'secondary_one_decay', 'experiment']='secondary', impact_parameter=(0., 3.), uball_multiplicity=(1, 25), histname=None, verbose=1):
        """ A class for handling TH2D tranverse momentum / A v.s. rapidity / beam-rapidity
        Parameters
        ----------
        particle : str
        path : str / pathlib.Path
            path of the root file containing the h2
        histname : str
            fullname of the h2 in `path`
        """
        self.particle = e15190.Particle(particle)
        self.reaction = reaction
        self.skyrme = 'EXP' if skyrme is None else skyrme
        self.impact_parameter = impact_parameter
        self.uball_multiplicity = uball_multiplicity
        self.mode = mode

        if path is None:
            if skyrme != 'EXP':
                path = pathlib.Path(self.DIR, f'{reaction}_{skyrme}.root')
            if not path.exists():
                raise IOError(f'file not found : {str(path)}')
            if not (self.path.is_file() and os.access(self.path, os.R_OK)):
                raise IOError(f'check file permission {str(path)}')

        if histname is None:
            if skyrme != 'EXP':
                histname = f'h2_pt_rapidity_{mode}_{particle}'
            else:
                histname = f'h2_pt_rapidity_{particle}'

        with root6.TFile(str(path)) as file:
            try:
                hist = file[histname]
            except:
                objname = file.keys()[0]
                if verbose == 1:
                    warnings.warn(
                        f'hist name {histname} not found in root file. We will use the first object in the file, i.e. `{objname}`')
                hist = file[objname]

            self.df = hist_reader.hist2d_to_df(hist)

    def query(self, df=None, xrange=(0.4, 0.6), yrange=(0., 600.), inplace=False):
        if df is None:
            df = self.df.copy()

        dx = np.diff(np.unique(df.x)[0:2])
        x1 = xrange[0] - dx/2.
        x2 = xrange[1] + dx/2
        dy = np.diff(np.unique(df.y)[0:2])
        y1 = yrange[0] - dy/2.
        y2 = yrange[1] + dy/2

        # In TH2D, projection is done on the bin range (bin1, bin2) instead of numerical range
        # e.g. if bin width is 0.005, then 0.396 is included in ROOT Projection.
        # 0.395 is not included becasue any value exactly the same as the interval would be filled in the bin on the left. so (> and <=)

        df.query(f'x > {x1} & x <= {x2} & y > {y1} & y <= {y2}', inplace=True)
        if inplace:
            self.df = df
        return df
    
    # def correct_coverage(self, df=None, xrange=(0.4, 0.6), yrange=(0., 600.)):
    #     df = self.query(df, xrange, yrange)
    #     df.query('z > 0.0', inplace=True)
    #     df.reset_index(inplace=True, drop=True)
    #     correction = np.ones(df.shape[0])
    #     for val, subdf in df.groupby('y'):
    #         if subdf.shape[0] == 1:
    #             correction[np.array(subdf.index)] = 0.
    #         else:
    #             dx = np.diff(np.array(subdf.x)[0:2])
    #             x1 = np.min(subdf.x) + dx * np.random.uniform(-0.5, 0.5)
    #             x2 = np.max(subdf.x) + dx * np.random.uniform(-0.5, 0.5)
    #             correction[np.array(subdf.index)] = abs(x2-x1)
    #     df['correction'] = correction

    #     df.query('correction > 0.0', inplace=True)
    #     df['correction'] = df['correction'] / np.diff(xrange)
    #     df['z'] = df['z']/df['correction']
    #     df['z_err'] = df['z_err']/df['correction']
    #     df.drop('correction', axis=1, inplace=True)
    #     return df

    def PtSpectrum(self, xrange=(0.4, 0.6), yrange=(0., 600.), bins=30,  correct_coverage=False, correct_range=(0., 600)):
        """ calculate ProjectionX of the h2
        Parameters
        ----------
        bins : int
            number of bins in Pt spectra
        correct_coverage : bool
            will correct for coverage if True
        """
        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=correct_range)
        else:
            df = self.df.query(f'x > {xrange[0]} & x <= {xrange[1]} & y > {yrange[0]} & y <= {yrange[1]}')

        df = df.query('z > 0.0')

        hist, edges = np.histogram(
            df.y, range=yrange, bins=bins, weights=df['z'])
        histerr, edges = np.histogram(
            df.y, range=yrange, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)

        norm = np.diff(yrange) / bins * np.diff(xrange)

        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist / norm,
            'y_err': histerr/norm,
            'y_ferr':  np.divide(histerr, hist, out=np.zeros_like(histerr), where=(hist != 0.0))
        })

    def RapidityDistribution(self, xrange=(0,1), yrange=(0,800), bins=100):
        """ Project to X-axis to get distribution of `rapidity_lab` / `beam_rapidity`
        """
        df = self.df.query(f'x > {xrange[0]} & x <= {xrange[1]} & y > {yrange[0]} & y <= {yrange[1]}')
        hist, _ = np.histogram(df['x'], weights=df['z'], bins=bins, range=xrange)
        hist_err, _ = np.histogram(df['x'], weights=df['z_err']**2, bins=bins, range=xrange)
        hist_err = np.sqrt(hist_err)
        
        x = np.linspace(*xrange, bins+1)
        return pd.DataFrame({
            'x' : 0.5 * (x[1:] + x[:-1]),
            'y' : hist,
            'y_err' : hist_err,
            'y_ferr' : np.divide(hist_err, hist, where=(hist!=0.))
        })

"""
    def RapidityCMS(self, range=(-0.5, 0.5), bins=100):
        df = self.df.copy()
        df.query('z > 0.0', inplace=True)
        dx = np.diff(np.unique(df.x)[0:2])
        x = df.x.to_numpy() + dx * np.random.uniform(-0.5, 0.5, size=len(df.x))

        rapidity_cms = x * self.beam_rapidity - 0.5 * \
            np.log((1+self.betacms)/(1-self.betacms))

        hist, edges = np.histogram(
            rapidity_cms, range=range, bins=bins, weights=df['z'])
        histerr, edges = np.histogram(
            rapidity_cms, range=range, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)

        norm = np.diff(range) / bins

        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist/norm,
            'y_err': histerr/norm,
            'y_ferr': np.divide(histerr, hist, out=np.zeros_like(histerr), where=(hist != 0.0))
        })

    def plotPtRapidity(self, ax=None, range=[[0., 1.], [0., 600.]], bins=[100, 600], **kwargs):
        df = self.df.copy()
        return df_helper.plot2d(ax, df, range=range, bins=bins, **kwargs)

    def plotPtSpectrum(self, ax=None, xrange=(0.4, 0.6), yrange=(0, 600), bins=30, correct_coverage=False, correct_range=(0., 600), **kwargs):
        df = self.PtSpectrum(xrange=xrange, yrange=yrange,
                             bins=bins, correct_coverage=correct_coverage, correct_range=correct_range)
        return df_helper.plot1d(ax, df, **kwargs)

    def plotRapidityCMS(self, ax=None, range=(-0.5, 0.5), bins=100, **kwargs):
        df = self.RapidityCMS(range=range, bins=bins)
        return df_helper.plot1d(ax, df, **kwargs)

    def EkinThetaCMS(self, xrange=(0.4, 0.6), yrange=(0., 600.), correct_coverage=False):
        df = self.correct_coverage(xrange=xrange, yrange=yrange) if correct_coverage else self.query(
            xrange=xrange, yrange=yrange)

        mass = 938.272

        dx = np.diff(np.unique(df.x)[0:2])
        dy = np.diff(np.unique(df.y)[0:2])

        pt = df.y + dy * np.random.uniform(-0.5, 0.5, size=len(df.y))
        rapidity_lab = df.x + dx * \
            (np.random.uniform(-0.5, 0.5, size=len(df.x)))
        rapidity_lab *= self.beam_rapidity
        rapidity_cms = rapidity_lab - 0.5 * \
            np.log((1.+self.betacms)/(1.-self.betacms))

        pzcms = np.sqrt(pt**2 + mass**2) * np.sinh(rapidity_cms)
        ekincms = np.sqrt(pt**2 + mass**2) * np.cosh(rapidity_cms) - mass
        thetacms = np.degrees(np.arctan2(pt, pzcms))

        return pd.DataFrame({
            'x': ekincms,
            'y': thetacms,
            'z': df['z'],
            'z_err': df['z_err'],
            'z_ferr': np.divide(df['z_err'], df['z'], out=np.zeros_like(df['z_err']), where=(df['z'] > 0.0))
        })

    def MomentumRapidityCMS(self, xrange=(0., 1.), yrange=(0., 600.), correct_coverage=False):
        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=yrange)
        else:
            df = self.query(xrange=xrange, yrange=yrange)

        mass = 938.272
        dx = np.diff(np.unique(df.x)[0:2])
        dy = np.diff(np.unique(df.y)[0:2])
        pt = df.y + dy * np.random.uniform(-0.5, 0.5, size=len(df.y))
        rapidity_lab = df.x + dx * \
            (np.random.uniform(-0.5, 0.5, size=len(df.x)))
        rapidity_lab *= self.beam_rapidity
        rapidity_cms = rapidity_lab - 0.5 * \
            np.log((1.+self.betacms)/(1.-self.betacms))

        pzcms = np.sqrt(pt**2 + mass**2) * np.sinh(rapidity_cms)
        pcms = np.sqrt(pt**2 + pzcms**2)

        return pd.DataFrame({
            'x': pcms,
            'y': rapidity_cms,
            'z': df['z'],
            'z_err': df['z_err'],
            'z_ferr': np.divide(df['z_err'], df['z'], out=np.zeros_like(df['z_err']), where=(df['z'] > 0.0))
        })

    def MomentumCMS(self, xrange=(0., 1.), yrange=(0., 800.), bins=40, cuts=[(-1.0, 1.0), (0., 800.)], correct_coverage=False):
        # Projecti Momentum from Momentum-Rapidity
        df = self.MomentumRapidityCMS(
            xrange=xrange, yrange=yrange, correct_coverage=correct_coverage)

        df.query(
            f'y >= {cuts[0][0]} & y <= {cuts[0][1]} & x >= {cuts[1][0]} & x <= {cuts[1][1]}', inplace=True)

        hist, edges = np.histogram(
            df.x, range=yrange, bins=bins, weights=df['z'])

        histerr, edges = np.histogram(
            df.x, range=yrange, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)
        norm = np.diff(yrange) / bins * np.diff(xrange)

        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist / norm,
            'y_err': histerr/norm,
            'y_ferr':  np.divide(histerr, hist, out=np.zeros_like(histerr), where=(hist != 0.0))
        })

    def Multiplicity(self, xrange=(0.4, 0.6), yrange=(0., 600.), bins=30,  correct_coverage=False, correct_range=(0., 600)):
        # return number of particle after normalization
        df = self.PtSpectrum(xrange=xrange, yrange=yrange, bins=bins,
                             correct_coverage=correct_coverage, correct_range=correct_range)
        return np.sum(df.y) * np.diff(xrange) * np.diff(yrange)/bins



class EkinThetaCMS:
    def __init__(self, particle, df=None, reaction=None, path=None, skyrme='SkM', impact_parameter=(0., 3.), uball_multiplicity=(1, 25), mode='all'):

        if df is None:
            df = PtRapidityLAB(
                particle=particle, reaction=reaction, path=path, skyrme=skyrme,
                impact_parameter=impact_parameter, uball_multiplicity=uball_multiplicity, mode=mode
            )
            df = df.EkinThetaCMS()
        df.columns = ['x', 'y', 'z', 'z_err', 'z_ferr']
        self.df = df

    def query(self, df=None, xrange=(0.4, 0.6), yrange=(0., 600.)):
        if df is None:
            df = self.df.copy()

        dx = np.diff(np.unique(df.x)[0:2])
        x1 = xrange[0] - dx/2.
        x2 = xrange[1] + dx/2
        dy = np.diff(np.unique(df.y)[0:2])
        y1 = yrange[0] - dy/2.
        y2 = yrange[1] + dy/2

        df.query(
            f'x >= {x1} & x <= {x2} & y >= {y1} & y <= {y2}', inplace=True)
        return df

    def EkinCMS(self, xrange=(0., 200.), yrange=(0., 180.), bins=50):
        df = self.query(xrange=xrange, yrange=yrange)
        df.query('z > 0.0', inplace=True)
        hist, edges = np.histogram(
            df.x, range=xrange, bins=bins, weights=df['z'])
        histerr, edges = np.histogram(
            df.x, range=xrange, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)

        norm = np.diff(xrange) / bins
        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist / norm,
            'y_err': histerr / norm,
            'y_ferr': np.divide(histerr, hist, out=np.zeros_like(histerr), where=(hist > 0.0))
        })

    def ThetaCMS(self, xrange=(0., 200.), yrange=(0., 180.), bins=180):
        df = self.query(xrange=xrange, yrange=yrange)
        df.query('z > 0.0', inplace=True)
        hist, edges = np.histogram(
            df.y, range=yrange, bins=bins, weights=df['z'])
        histerr, edges = np.histogram(
            df.y, range=yrange, bins=bins, weights=df['z_err']**2)
        histerr = np.sqrt(histerr)

        norm = np.diff(xrange) / bins
        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist / norm,
            'y_err': histerr / norm,
            'y_ferr': np.divide(histerr, hist, out=np.zeros_like(histerr), where=(hist > 0.0))
        })

"""