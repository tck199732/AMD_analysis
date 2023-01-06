import pathlib
import os
import iminuit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

from pyamd import PROJECT_DIR
from pyamd.e15190 import e15190
from pyamd.utilities import root6, dataframe


hist_reader = root6.HistogramReader()
df_helper = dataframe.DataFrameHelper()

class EmissionTimeFile:
    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()

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


class EmissionTime:
    DIR = f'{str(PROJECT_DIR)}/result/emission_time'

    def __init__(self, particle, path=None, reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(0., 3.)):
        self.particle = e15190.particle(particle)
        if path is None:
            path = f'{self.DIR}/{reaction}_{skyrme}_bmin{impact_parameter[0]:.1f}_bmax{impact_parameter[1]:.1f}.root'

        spectrafile = EmissionTimeFile(path)
        self.df = spectrafile.get_histogram(
            keyword=f'h2_time_momentum_{particle}')
        self.profile = self.df.ProfileX(
            name='profile', firstybin=1, lastybin=-1)
        self.profile = hist_reader.hist1d_to_df(self.profile)
        self.df = hist_reader.hist2d_to_df(self.df)

    '''
    so far option g and s are correct, default off by a factor of \sqrt(# of events))
    
    def AverageEmissionTime(self, range=(0, 600), bins=30, drop_zeros=True, option=None):
        df = self.df.copy()
        if drop_zeros:
            df.query('z > 0.0', inplace=True)
        df.query(
            f'y >= {range[0]} & y <= {range[1]}', inplace=True)

        hist = df.groupby(['x']).apply(
            lambda x: np.average(x.y, weights=x['z']))
        
        if option is None:
            stderr= lambda x: np.sqrt((np.average(x['y']**2, weights=x['z']) - np.average(x['y'], weights=x['z'])**2)/np.sum(x['z'])**2 * np.sum(x['z]**2))
            hist_err = df.groupby(['x']).apply(stderr)
        elif option == 's':
            stderr= lambda x: np.sqrt((np.average(x['y']**2, weights=x['z']) - np.average(x['y'], weights=x['z'])**2))
            hist_err = df.groupby(['x']).apply(stderr)
        elif option == 'g':
            hist_err = df.groupby(['x']).apply(lambda x: 1./np.sqrt(np.sum(x['z_err'])))
        
        x = hist.index.to_numpy()
        hist = hist.to_numpy()
        hist_err = hist_err.to_numpy()

        hist, x_edges = np.histogram(
            x, range=range, bins=bins, weights=hist)
        hist_err, x_edges = np.histogram(
            x, range=range, bins=bins, weights=hist_err**2)
        hist_err = np.sqrt(hist_err)
        
        dt = np.diff(range) / bins
        sumw, x_edges = np.histogram(x, range=range, bins=bins, weights=sumw)

        return pd.DataFrame({
            'x': 0.5 * (x_edges[1:] + x_edges[:-1]),
            'y': hist / dt,
            'y_err': hist_err / dt,
            'y_ferr': np.divide(hist_err, hist, where=(hist > 0.0), out=np.zeros_like(hist_err))
        })
    '''

    def AverageEmissionTime(self, range=(0, 600), bins=30, drop_zero_time=True):
        df = self.profile.copy()
        if drop_zero_time:
            df.query('y >= 1.', inplace=True)
        hist, x_edges = np.histogram(
            df.x, range=range, bins=bins, weights=df.y)
        hist_err, x_edges = np.histogram(
            df.x, range=range, bins=bins, weights=df['y_err']**2)
        hist_err = np.sqrt(hist_err)
        dt = np.diff(range) / bins
        return pd.DataFrame({
            'x': 0.5 * (x_edges[1:] + x_edges[:-1]),
            'y': hist / dt,
            'y_err': hist_err / dt,
            'y_ferr': np.divide(hist_err, hist, where=(hist > 0.0), out=np.zeros_like(hist_err))
        })

    def plot2d(self, ax=None, df=None, cmap='jet',  drop_zero_time=True, **kwargs):
        cmap = copy(plt.cm.get_cmap(cmap))
        cmap.set_under('white')

        kw = dict(
            cmap=cmap,
            range=[(0, 600), (0, 500)],
            bins=[30, 100]
        )
        kw.update(kwargs)

        if df is None:
            df = self.df.copy()

        if drop_zero_time:
            df.query('y >= 1.0', inplace=True)

        if ax is None:
            ax = plt.gca()
        ax.hist2d(df.x, df.y, weights=df['z'], **kw)
        return ax

    def plot1d(self, ax=None, range=(0, 600), bins=30, drop_large_error=True, thres=0.05, **kwargs):
        kw = dict(
            fmt='.'
        )
        kw.update(kwargs)
        df = self.AverageEmissionTime(range=range, bins=bins)
        df.query('y > 0.0', inplace=True)
        # if drop_large_error:
        # df.query(f'y_ferr < {thres}', inplace=True)
        if ax is None:
            ax = plt.gca()

        ax.errorbar(df.x, df.y, yerr=df['y_err'], **kw)
        return ax


class EmissionPosition:
    def __init__(self, path, particle, reaction, skyrme):
        self.reaction = e15190.reaction(reaction)
        self.particle = e15190.particle(particle)
        self.skyrme = skyrme

        spectrafile = EmissionTimeFile(path)
        df = spectrafile.get_histogram(keyword=f'h2_time_position_{particle}')
        # x = position; y = time;
        self.df = hist_reader.hist2d_to_df(df)

    def get_source_fcn(self, bins=400, range=(0, 40), cuts=[(0, 40), (1., 500)], normalize=True, jacobian=True):

        df = self.df.copy()
        df.query(
            f'x >= {cuts[0][0]} & x <= {cuts[0][1]} & y >= {cuts[1][0]} & y<= {cuts[1][1]}', inplace=True)

        x = np.linspace(*range, bins+1)
        x = 0.5 * (x[1:] + x[:-1])

        df['weight'] = df['z']
        df['err_weight'] = df['z_err']
        if jacobian:
            df['weight'] = df['z'] / (4. * np.pi * df['x']**2)
            df['err_weight'] = df['z_err'] / (4. * np.pi * df['x']**2)

        hist, edges = np.histogram(
            df.x, bins=bins, range=range, weights=df['weight'])
        histerr, edges = np.histogram(
            df.x, bins=bins, range=range, weights=df['err_weight'])
        histerr = np.sqrt(histerr)

        if normalize:
            scale = np.sum(4. * np.pi * x**2 * hist * (x[1] - x[0]))
            hist /= scale
            histerr /= scale

        return pd.DataFrame({
            'x': x,
            'y': hist,
            'y_err': histerr,
            'y_ferr': np.divide(histerr, hist, where=hist != 0., out=np.zeros_like(histerr)),
        })


    def get_source_size(self, bins=400, range=(0, 40), cuts=[(0, 40), (0, 500)], ):
        df = self.get_source_fcn(
            bins=bins, range=range, cuts=cuts, normalize=True, jacobian=True)
        half_maximum = df.y.max() / 2.
        id = np.abs(df.y - half_maximum).arcmin()
        return df.x[id]

    # should only run once.
    def fit_gaussian(self, bins=400, range=(0, 40), cuts=[(0, 40), (0, 500)]):
        df = self.get_source_fcn(
            bins=bins, range=range, cuts=cuts, normalize=True, jacobian=True)
        
        global model
        model = lambda r, l, r0 : l / (2 * np.sqrt(np.pi) * r0)**3 * np.exp(-r**2 / 4 / r0**2)

        least_squares = iminuit.cost.LeastSquares(df.x, df.y, df['y_err'], model)
        minuit = iminuit.Minuit(least_squares, l=1.0, r0=3.)
        minuit.migrad()
        return minuit
        # return (*minuit.values, model(df.x, *minuit.values))

    def get_gaussian_parameters(self):
        return self.fit_gaussian().values
    def get_gaussian_errors(self):
        return self.fit_gaussian().errors
    
        
    
    def plot2d(self, ax=None, df=None, cmap='jet',  drop_zero_time=True, **kwargs):
        cmap = copy(plt.cm.get_cmap(cmap))
        cmap.set_under('white')

        kw = dict(
            cmap=cmap,
            range=[(0, 40), (0, 500)],
            bins=[400, 100]
        )
        kw.update(kwargs)

        if df is None:
            df = self.df.copy()

        if drop_zero_time:
            df.query('y >= 1.0', inplace=True)

        if ax is None:
            ax = plt.gca()
        ax.hist2d(df.x, df.y, weights=df['z'], **kw)
        return ax

    def plot1d(self, ax=None, bins=400, range=(0, 40), cuts=[(0, 40), (0, 500)], normalize=True, jacobian=True, **kwargs):
        df = self.get_source_fcn(bins=bins, range=range, cuts=cuts, normalize=normalize, jacobian=jacobian)
        return df_helper.plot1d(ax, df, **kwargs)
        
        