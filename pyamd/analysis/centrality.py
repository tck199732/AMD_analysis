# plotting
import matplotlib.pyplot as plt

# I/O
import pathlib
import warnings

# analysis
import numpy as np
import pandas as pd
from astropy import units
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.optimize import curve_fit

# miscellaneous
from copy import copy
from functools import cache

# local
from pyamd import PROJECT_DIR
from pyamd.utilities import dataframe, minuit, root6, symwrap

hist_reader = root6.HistogramReader()
df_helper = dataframe.DataFrameHelper()


class Multiplicity_ImpactParameter:
    DIR = f'{PROJECT_DIR}/result/centrality/b10fm'

    def __init__(self, path=None, histname='h2_multi_b', reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(0., 10.), mode='3', verbose=0):
        """ A Class handle calculation from 2D histogram Impact-Parameter VS Multiplicity
        Parameter
        ---------
        path : str, pathlib.Path
            if None, use the default path
        histname : str
            name of the TH2D in the root file
        reaction : str
            collision string in the format of {beam}{target}E{energy}
        skyrme : str
            skyrme parameter set used in the simulation. For now, either `SkM` or `SLy4` or `SLy4_L108`
        impact_parameter : tuple of float
            range of b in the simulation. This is only included for completeness.
        mode : str, int, float
            simulation mode : `21` means data from table21. `3` means data from table3 (sec. decay)
        verbose : int 
            `0` = quiet
            `1` = issue warning in case histogram is not found
        """
        self.reaction = reaction
        self.skyrme = skyrme
        self.impact_parameter = impact_parameter
        self.mode = str(mode)

        if path is None:
            path = pathlib.Path(
                f'{self.DIR}/{reaction}_{skyrme}_table{mode}.root')
            if not path.exists():
                raise ValueError(f'path does not exist : {str(path)}')

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

    def MultiplicitySpectra(self, range=(-0.5, 24.5), cut=(0., 10.), bins=25, normalize=True):
        """ Projection to get multiplicity spectra
        Parameter
        ---------
        range : tuple of float
            range of multiplicity
        cut : tuple of float
            cut on impact-parameter
        bins : int
            number of bins
        normalize : bool
            if True, normalize by bin width
        """
        df = self.df.copy()
        binwidth = np.diff(df['y'][0:2])
        cut1 = cut[0] - binwidth/2.
        cut2 = cut[1] + binwidth/2.

        df.query('z > 0 ', inplace=True)
        df.query(f'y >= {float(cut1)} & y < {float(cut2)}', inplace=True)

        hist, bin_edges = np.histogram(
            df['x'], weights=df['z'], range=range, bins=bins)
        hist_err, bin_edges = np.histogram(
            df['x'], weights=df['z_err']**2, range=range, bins=bins)
        hist_err = np.sqrt(hist_err)

        # normalization
        dM = np.diff(range)/bins
        norm = np.abs(dM) if normalize else 1.

        return pd.DataFrame({
            'x': 0.5 * (bin_edges[1:] + bin_edges[:-1]),
            'y': hist / norm,
            'y_err': hist_err / norm,
            'y_ferr': np.divide(hist_err, hist, out=np.zeros_like(hist_err), where=(hist > 0.0))
        })

    def ImpactParameterSpectra(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, normalize=False):
        """ Projection to get impact-parameter spectra
        Parameter
        ---------
        range : tuple of float 
            range of b
        cut : tuple of float 
            cut on multiplity
        bins : int
            number of bins
        normalize : bool
            if True, normalize by bin width
        """
        df = self.df.copy()
        df.query(f'x >= {cut[0]} & x <= {cut[1]}', inplace=True)

        hist, bin_edges = np.histogram(
            df['y'], weights=df['z'], range=range, bins=bins)
        hist_err, bin_edges = np.histogram(
            df['y'], weights=df['z_err']**2, range=range, bins=bins)
        hist_err = np.sqrt(hist_err)

        # normalization
        db = np.diff(range)/bins
        norm = np.abs(db) if normalize else 1.

        return pd.DataFrame({
            'x': 0.5 * (bin_edges[1:] + bin_edges[:-1]),
            'y': hist / norm,
            'y_err': hist_err / norm,
            'y_ferr': np.divide(hist_err, hist, out=np.zeros_like(hist_err), where=(hist > 0.0))
        })

    def DifferentialCrossSection(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=None):
        """ Obtain differential cross section :math: `d\sigma/db`. Note that we assume the input histogram is normalized by **`pure simulation events`** not **`filtered events`**.
        Parameter
        ---------
        range : tuple of float
            range of b
        cut : tuple of float
            cut on multiplicity
        bins : int
            number of bins
        unit : str
            `None` : use the same unit as in simulation, i.e. fm^2
            `b` : barn
            `mb` : 0.1 barn
        """
        # normalized impact parameter spectra dM/db
        df = self.ImpactParameterSpectra(
            range=range, cut=cut, bins=bins, normalize=True)

        # first get total cross section \pi bmax^2 (fm^2)
        cs = np.pi * range[1]**2 * units.fm**2
        if unit is None:
            pass
        elif unit == 'mb':
            cs = cs.to(1e-3 * units.barn)
        elif unit == 'b':
            cs = cs.to(units.barn)

        # dsigma/db = \pi bmax^2 * normalized dM/db
        dsigma_db = cs.value * df.y
        error = cs.value * df['y_err']

        df = pd.DataFrame({
            'x': df.x,
            'y': dsigma_db,
            'y_err': error
        })
        df_helper.calculate_relative_error(df)
        return df

    def CrossSection(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=None):
        """ Integrate :math: `d\sigma/db` to get :math: `\sigma`
        range : tuple of float
            range of b
        cut : tuple of float
            cut on multiplicity
        bins : int
            number of bins
        unit : str
            `None` : use the same unit as in simulation, i.e. fm^2
            `b` : barn
            `mb` : 0.1 barn
        """
        df = self.DifferentialCrossSection(
            range, cut, bins, unit)  # full range
        db = np.diff(df.x[0:2])
        return np.sum(df['y']) * db, np.sqrt(np.sum(df['y_err'])) * db

    def BmaxFromCrossSection(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=None):
        sig, err = self.CrossSection(range, cut, bins, unit)
        return np.sqrt(sig / np.pi), err / 2. / np.sqrt(np.pi * sig)

    def Mapping(self, ranges=[[-0.5, 24.5], [0., 10.]], bins=[25, 20], unit=None):
        """ Calculate the mapping from charged-particle multiplicty to Impact-Parameter
        Parameter
        ---------
        range: 
            range of multiplicty (0) and b (1)
        bins:
            bins of multiplicty (0) and b (1) 
        unit:
            unit of impact-parameter
            `None` : fm
        """
        df = self.MultiplicitySpectra(
            range=ranges[0], cut=ranges[1], bins=bins[0])

        # in case the spectra isn't normalized
        df['y'] /= np.sum(df['y'])
        df['y_err'] /= np.sum(df['y'])

        bmax, bmax_err = self.BmaxFromCrossSection(
            range=ranges[1], cut=ranges[0], bins=bins[1], unit=unit)

        y = bmax * np.sqrt(1. - np.cumsum(df['y']))
        yerr = np.sqrt(
            (0.5 * bmax**2 / y * np.sqrt(np.cumsum(df['y_err']**2))) ** 2. +
            (bmax_err * y / bmax) ** 2.
        )

        df = pd.DataFrame({
            'x': df.x.to_numpy(),
            'y': y,
            'y_err': yerr,
        })

        df_helper.calculate_relative_error(df, inplace=True)
        return df

    @staticmethod
    def create_model_dsigma_db():
        expr = 'a * x / (1 + exp((x - b) / c))'
        vars = ['x', 'a', 'b', 'c']
        return symwrap.expression.make_function(vars, expr)

    @staticmethod
    def create_model_dsigma_db_error():
        expr = 'a * x / (1 + exp((x - b) / c))'
        vars = ['x', 'a', 'b', 'c']
        return symwrap.expression.error_propagation(vars, pars=['x'], expr=expr)

    @staticmethod
    def model_dsigma_db(x, norm, b0, db):
        f = Multiplicity_ImpactParameter.create_model_dsigma_db()
        return f(x, norm, b0, db)

    @staticmethod
    def model_dsigma_db_error(x, norm, b0, db, norm_err, b0_err, db_err):
        ferr = Multiplicity_ImpactParameter.create_model_dsigma_db_error()
        return ferr(x, norm, b0, db, norm_err, b0_err, db_err)

    def fit_dsigma_db_iminuit(self, df=None, range=(0., 10.), bins=20, unit=None, verbose=1, output_info=False):
        if df is None:
            df = self.DifferentialCrossSection(
                range=range, bins=bins, unit=None)

        cost_fcn = LeastSquares(
            df.x, df.y, df['y_err'], Multiplicity_ImpactParameter.model_dsigma_db)
        m = Minuit(cost_fcn, norm=2.*np.pi, b0=9., db=0.1)
        m.migrad()
        m.hesse()

        fit_info = {
            'chi2': m.fval,
            'ndof': df.x.shape[0] - m.nfit,
        }
        fit_info.update({
            p : v for p, v in zip(m.parameters, m.values)
        })

        fit_info.update({
            f'{p}_err' : e for p, e in zip(m.parameters, m.errors)
        })

        if verbose != 0:
            print('\n'.join(
                [f'chi^2 / dof = {m.fval:.1f} / {len(df.x) - m.nfit}'] +
                [f'{p} = {v:.3f} +/- {e:.3f}' for p, v,
                 e in zip(m.parameters, m.values, m.errors)]
            ))

        x = df.x.to_numpy()
        y = Multiplicity_ImpactParameter.model_dsigma_db(
            x, *m.values) * units.fm**2
        yerr = Multiplicity_ImpactParameter.model_dsigma_db_error(
            x, *m.values, *m.errors) * units.fm**2

        if unit == 'mb':
            y = y.to(1e-3 * units.barn)
            yerr = yerr.to(1e-3 * units.barn)

        elif unit == 'b':
            y = y.to(units.barn)
            yerr = yerr.to(units.barn)

        df = pd.DataFrame({
            'x': x,
            'y': y.value,
            'y_err': yerr.value
        })
        df_helper.calculate_relative_error(df)

        if output_info:
            return df, fit_info
        return df

    def fit_dsigma_db_TMinuit(self, df=None, range=(0., 10.), bins=20, unit=None):
        """ Fit the differential cross-section using TMinuit.
        Parameter
        ---------
        df : pd.DataFrame
            dataframe of the :math: `d\sigma/db`
        range : tuple of float
            range of b used in calculating :math: `d\sigma/db`
        bins : int
            number of bins used in calculating :math: `d\sigma/db`
        unit : str
            This unit is used in the final output. The fit is always done using the origianl unit of :math: `d\sigma/db`, i.e. :math: `fm^2` or `unit = None`
            `mb` : 1e-1 barn
        """
        if df is None:
            df = self.DifferentialCrossSection(
                range=range, bins=bins, unit=None)
        global fcn
        global xarray
        xarray = df.x.to_numpy()

        def fcn(npar, gin, f, par, iflag):
            chisq = (Multiplicity_ImpactParameter.model_dsigma_db(xarray, par[0], par[1], par[2]) -
                     df.y.to_numpy()) / df['y_err']
            chisq = np.sum(chisq**2)
            f.value = chisq

        gMinuit = minuit.TMinuit(3, fcn)
        gMinuit.set_verbose(0)
        gMinuit.set_parameter(0, 'norm', 2. * np.pi, 0.1)
        gMinuit.set_parameter(1, 'b0', 9., 0.1)
        gMinuit.set_parameter(2, 'db', 0.1, 0.01)
        gMinuit.fit()

        norm = gMinuit.get_parameter('norm')
        b0 = gMinuit.get_parameter('b0')
        db = gMinuit.get_parameter('db')

        norm_err = gMinuit.get_error('norm')
        b0_err = gMinuit.get_error('b0')
        db_err = gMinuit.get_error('db')

        x = df.x.to_numpy()
        y = Multiplicity_ImpactParameter.model_dsigma_db(
            x, norm, b0, db
        ) * units.fm**2
        yerr = Multiplicity_ImpactParameter.model_dsigma_db_error(
            x, norm, b0, db, norm_err, b0_err, db_err
        ) * units.fm**2

        if unit == 'mb':
            y = y.to(1e-3 * units.barn)
            yerr = yerr.to(1e-3 * units.barn)

        elif unit == 'b':
            y = y.to(units.barn)
            yerr = yerr.to(units.barn)

        return pd.DataFrame({
            'x': x,
            'y': y.value,
            'y_err': yerr.value,
            'y_ferr': 0.
        })

    def plot_multiplicity(self, ax=None, range=(-0.5, 24.5), cut=(0., 10.), bins=25, **kwargs):
        df = self.MultiplicitySpectra(range=range, bins=bins, cut=cut)
        kw = dict(fmt='.')
        kw.update(kwargs)
        return df_helper.plot1d(ax, df, **kw)

    def plot_impact_parameter(self, ax=None, range=(0., 10.), cut=(-0.5, 24.5), bins=20, **kwargs):
        df = self.ImpactParameterSpectra(range=range, bins=bins, cut=cut)
        kw = dict(
            fmt='.',
        )
        kw.update(kwargs)
        return df_helper.plot1d(ax, df, **kw)

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

        df = df.query('y > 0. & y_err > 0.')

        par, cov = curve_fit(
            # lambda x, *par : model(x, *par),
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
        return df_helper.plot1d(ax, df, **kw)

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
