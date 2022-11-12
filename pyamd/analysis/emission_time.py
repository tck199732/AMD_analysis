import pathlib
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
from pyamd import PROJECT_DIR
from pyamd.utilities import root6
hist_reader = root6.HistogramReader()


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
        self.particle = particle
        if path is None:
            path = f'{self.DIR}/{reaction}_{skyrme}_bmin{impact_parameter[0]:.1f}_bmax{impact_parameter[1]:.1f}.root'

        spectrafile = EmissionTimeFile(path)
        self.df = spectrafile.get_histogram(keyword=f'h2_time_momentum_{particle}')
        self.profile = self.df.ProfileX(name='profile', firstybin=1, lastybin=-1)
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

    def AverageEmissionTime(self, range=(0, 600), bins=30):
        df = self.profile.copy()
        hist, x_edges = np.histogram(df.x, range=range, bins=bins, weights=df.y)
        hist_err, x_edges = np.histogram(df.x, range=range, bins=bins, weights=df['y_err']**2)
        hist_err = np.sqrt(hist_err)
        dt = np.diff(range) / bins
        return pd.DataFrame({
            'x': 0.5 * (x_edges[1:] + x_edges[:-1]),
            'y': hist / dt,
            'y_err': hist_err / dt,
            'y_ferr': np.divide(hist_err, hist, where=(hist > 0.0), out=np.zeros_like(hist_err))
        })


    def plot2d(self, ax=None, df=None, cmap='jet', **kwargs):
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

        if ax is None:
            ax = plt.gca()
        ax.hist2d(df.x, df.y, weights=df['z'], **kw)
        return ax

    def plot1d(self, ax=None, range=(0, 600), bins=30, drop_large_error=True, thres=0.05, **kwargs):
        kw = dict(
            fmt = '.'
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