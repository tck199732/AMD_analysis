import pathlib
import os
import sys
import functools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy import modeling
import subprocess
from copy import copy

import pyamd
from pyamd import PROJECT_DIR
from pyamd.utilities import root6, helper, minuit
from pyamd.e15190 import e15190

hist_reader = root6.HistogramReader()
helper = helper.helper()


class MultiplicityFile:
    def __init__(self, path):
        self.path = pathlib.Path(__file__).parent.resolve()
        self.path = pathlib.Path(self.path, path).resolve()

    def get_names(self):
        if not self.path.exists():
            raise OSError(f'file not found : {str(self.path)}')
        if self.path.is_file() and os.access(self.path, os.R_OK):
            with root6.TFile(str(self.path)) as file:
                return [key.GetName() for key in file.GetListOfKeys()]

    def get_histogram(self, name=None, keyword=None):
        if name is None:
            names = self.get_names()
            names = [name for name in names if keyword.lower() in name.lower()]
            print(names)
            if len(names) > 1:
                raise ValueError('more than 1 objects is returned.')
            name = max(names, key=len)

        with root6.TFile(self.path) as file:
            hist = file.Get(name)
            hist.SetDirectory(0)  # will not be deleted when file is closed
            return hist


class Multiplicity_ImpactParameter:
    DIR = f'{PROJECT_DIR}/multiplicity_analysis/result'

    def __init__(self, path=None, reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(0., 10.), mode='all'):

        self.betacms = e15190.reaction(reaction).get_betacms()
        self.beam_rapidity = e15190.reaction(reaction).get_rapidity_beam()
        self.reaction = reaction

        if path is None:
            path = f'{self.DIR}/amd_{reaction}_{skyrme}_bmin{impact_parameter[0]:.1f}_bmax{impact_parameter[1]:.1f}.root'

        spectrafile = MultiplicityFile(path)
        self.df = spectrafile.get_histogram(keyword=f'h2_multi_b_{mode}')
        self.df = hist_reader.hist2d_to_df(self.df)

    def MultiplicitySpectra(self, range=(-0.5, 24.5), cut=(0., 10.), bins=25):

        df = self.df.copy()
        binwidth = np.diff(df['y'][0:2])
        cut1 = cut[0] - binwidth/2.
        cut2 = cut[1] + binwidth/2.
        df.query(f'y >= {cut1} & y < {cut2}', inplace=True)

        hist, bin_edges = np.histogram(
            df['x'], weights=df['z'], range=range, bins=bins)
        hist_err, bin_edges = np.histogram(
            df['x'], weights=df['z_err']**2, range=range, bins=bins)
        hist_err = np.sqrt(hist_err)

        # normalization
        dM = np.diff(range)/bins
        norm = np.abs(dM)
        return pd.DataFrame({
            'x': 0.5 * (bin_edges[1:] + bin_edges[:-1]),
            'y': hist / norm,
            'y_err': hist_err / norm,
            'y_ferr': np.divide(hist_err, hist, out=np.zeros_like(hist_err), where=(hist > 0.0))
        })

    def ImpactParameterSpectra(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20):

        df = self.df.copy()
        df.query(f'x >= {cut[0]} & x <= {cut[1]}', inplace=True)

        hist, bin_edges = np.histogram(
            df['y'], weights=df['z'], range=range, bins=bins)
        hist_err, bin_edges = np.histogram(
            df['y'], weights=df['z_err']**2, range=range, bins=bins)
        hist_err = np.sqrt(hist_err)

        # normalization
        db = np.diff(range)/bins
        norm = np.abs(db)

        return pd.DataFrame({
            'x': 0.5 * (bin_edges[1:] + bin_edges[:-1]),
            'y': hist / norm,
            'y_err': hist_err / norm,
            'y_ferr': np.divide(hist_err, hist, out=np.zeros_like(hist_err), where=(hist > 0.0))
        })

    def DifferentialCrossSection(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=1.):
        # normalized dM/db
        df = self.ImpactParameterSpectra(
            range=range, cut=cut, bins=bins)

        # dsigma/db = \pi bmax^2 * normalized dM/db
        cs = np.pi * range[1]**2
        if unit == 'mb':
            unit = 10.

        return pd.DataFrame({
            'x': df.x,
            'y': unit * cs * df.y,
            'y_err': unit * cs * df['y_err'],
            'y_ferr': df['y_ferr']
        })

    def CrossSection(self, unit=1.):
        df = self.DifferentialCrossSection()  # full range
        db = np.diff(df.x[0:2])
        if unit == 'mb':
            unit = 1./10
        elif unit == 'b':
            unit = 1./100

        return np.sum(df['y']) * db * unit

    def BmaxFromCrossSection(self, unit=1.):
        return np.sqrt(self.CrossSection(unit=unit) / np.pi)

    def Mapping(self, unit=1.):

        df = self.MultiplicitySpectra()
        df.query('y > 0.0', inplace=True)
        df['y'] = df['y'] / np.sum(df['y'])
        df['y_err'] = df['y_err'] / np.sum(df['y'])
        bmax = self.BmaxFromCrossSection(unit=unit)
        y = bmax * np.sqrt(np.abs(1. - np.cumsum(df['y'])))
        yerr = 0.5 * bmax**2 / y * np.sqrt(np.cumsum(df['y_err']**2))

        return pd.DataFrame({
            'x': df.x.to_numpy(),
            'y': y,
            'y_err': 0.,
            'y_ferr': 0.
            # 'y_err': yerr,
            # 'y_ferr': np.divide(yerr, y, out=np.zeros_like(yerr), where=(y > 0.0))
        })

    def plot_multiplicity(self, ax=None, range=(-0.5, 24.5), cut=(0., 10.), bins=25, **kwargs):
        df = self.MultiplicitySpectra(range=range, bins=bins, cut=cut)
        kw = dict(fmt='.')
        kw.update(kwargs)
        return helper.plot1d(ax, df, **kw)

    def plot_impact_parameter(self, ax=None, range=(0., 10.), cut=(-0.5, 24.5), bins=20, **kwargs):
        df = self.ImpactParameterSpectra(range=range, bins=bins, cut=cut)
        kw = dict(
            fmt='.',
        )
        kw.update(kwargs)
        return helper.plot1d(ax, df, **kw)

    def model_dsigma_db(self, x, norm, b0, db):
        return norm * x / (1 + np.exp((x-b0)/db))

    def fit_dsigma_db(self, df=None, range=(0., 10.), bins=20):
        if df is None:
            df = self.DifferentialCrossSection(
                range=range, bins=bins, unit=1.)
        global fcn
        global xarray
        xarray = df.x.to_numpy()

        def fcn(npar, gin, f, par, iflag):
            chisq = (self.model_dsigma_db(xarray, par[0], par[1], par[2]) -
                     df.y.to_numpy()) / df['y_err']
            chisq = np.sum(chisq**2)
            f.value = chisq

        gMinuit = minuit.TMinuit(3, fcn)
        gMinuit.set_parameter(0, 'norm', 6.28, 0.1)
        gMinuit.set_parameter(1, 'b0', 9., 0.1)
        gMinuit.set_parameter(2, 'db', 0.1, 0.01)
        gMinuit.fit()

        norm = gMinuit.get_parameter('norm')
        b0 = gMinuit.get_parameter('b0')
        db = gMinuit.get_parameter('db')

        return pd.DataFrame({
            'x': df.x.to_numpy(),
            'y': self.model_dsigma_db(df.x.to_numpy(), norm, b0, db),
            'y_err': 0.,
            'y_ferr': 0.
        })

    def plot_fitted_dsigma_db(self, ax=None, df=None, range=(0., 10.), bins=20, unit=1., **kwargs):
        df = self.fit_dsigma_db(df=df, range=range, bins=bins)

        unit = 10. if unit == 'mb' else 1.
        df['y'] = df['y'] * unit
        df['y_err'] = df['y_err'] * unit
        kw = dict(fmt='--')
        kw.update(kwargs)
        return helper.plot1d(ax, df, **kw)

    def plot_dsigma_db(self, ax=None, df=None, range=(0., 10.), bins=20, unit=1., **kwargs):
        if df is None:
            df = self.DifferentialCrossSection(
                range=range, bins=bins, unit=unit)
        kw = dict(fmt='.',)
        kw.update(kwargs)
        return helper.plot1d(ax, df, **kw)

    def plot2d(self, ax=None, df=None, cmap='jet', drop_zeros=True, **kwargs):
        if df is None:
            df = self.df.copy()

        if drop_zeros:
            df = df.query(f'z > 0.0')

        if ax is None:
            ax = plt.gca()

        cmap = copy(plt.cm.get_cmap(cmap))
        cmap.set_under('white')

        kw = dict(
            cmap=cmap,
            range=[[-0.5, 24.5], [0., 10.]],
            bins=[25, 100]
        )
        kw.update(kwargs)

        return ax.hist2d(df['x'], df['y'], weights=df['z'], **kw)

    def fit1d(self, df):

        def model(x, norm, mean, sigma):
            return norm * np.exp(-(x-mean)**2 / 2. / sigma**2)

        df.query('y > 0.0 & y_err > 0.0', inplace=True)

        par, cov = curve_fit(
            #lambda x, *par : model(x, *par),
            model,
            df.x.to_numpy(),
            df.y.to_numpy(),
            sigma=df['y_err'].to_numpy(),
            absolute_sigma=True
        )
        err = np.sqrt(np.diag(cov))
        return par, err

    def plotfit(self, ax=None, df=None, range=(-0.5, 24.5), bins=25., **kwargs):

        par, err = self.fit1d(df)
        x = np.linspace(*range, bins+1)
        norm, mean, sigma = par
        norm_err, mean_err, sigma_err = err
        y = norm * np.exp(-(x-mean)**2 / 2. / sigma**2)
        yerr = 0.
        kw = dict(fmt='--')
        kw.update(kwargs)

        df = pd.DataFrame({
            'x': x,
            'y': y,
            'y_err': yerr,
            'y_ferr': 0.
            # 'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })
        return helper.plot1d(ax, df, **kw)

    def plot_fitted_multiplicity(self, ax=None, df=None, cut=(0., 10.), range=(-0.5, 24.5), bins=25, **kwargs):

        if df is None:
            df = self.MultiplicitySpectra(cut=cut, range=range, bins=bins)
        kw = dict(fmt='--')
        kw.update(kwargs)
        return self.plotfit(ax, df, range=range, bins=bins, **kw)

    def plot_fitted_impact_parameter(self, ax=None, df=None, cut=(1, 25), range=(0., 10.), bins=20, **kwargs):
        if df is None:
            df = self.ImpactParameterSpectra(
                cut=cut, range=range, bins=bins)
        kw = dict(fmt='--')
        kw.update(kwargs)
        return self.plotfit(ax, df, range=range, bins=bins, **kw)

    def plot_mapping(self, ax=None, unit=1., **kwargs):
        df = self.Mapping(unit=unit)

        kw = dict(fmt='ko')
        kw.update(kwargs)
        return helper.plot1d(ax, df, **kw)

    def output(self, df, path, **kwargs):
        kw = dict(
            sep='\t',
            columns=df.columns,
            mode='w',
        )
        kw.update(kwargs)
        df.to_csv(path, kw)
