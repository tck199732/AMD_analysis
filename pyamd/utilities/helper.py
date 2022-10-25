import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import os
from copy import copy
from pyamd.utilities import root6


class helper:
    def __init__(self):
        pass
    # df2 / df1

    def ratio1d(self, df1, df2):
        df1 = df1.copy()
        df2 = df2.copy()
        if np.any(df1.x.to_numpy() != df2.x.to_numpy()):
            raise ValueError('dfs have differnet x.')
        x = df1.x
        y = np.where((df1.y != 0.0), df2.y / df1.y, 0.0)
        yerr = y * np.sqrt(df1['y_ferr']**2 + df2['y_ferr']**2)

        return pd.DataFrame({
            'x': x,
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })

    def plot1d(self, ax=None, df=None, drop_zeros=True, drop_large_err=False, rel_err=0.05, **kwargs):
        if ax is None:
            ax = plt.gca()
        if df is None:
            raise ValueError('Empty df')
        df = df.copy()
        if drop_zeros:
            df.query('y != 0.0', inplace=True)
        if drop_large_err:
            df.query(f'y_ferr < {rel_err}', inplace=True)

        kw = dict(
            fmt='.'
        )
        kw.update(kwargs)
        ax.errorbar(df.x, df.y, yerr=df['y_err'], **kw)
        return ax

    def plot2d(self, ax=None, df=None, cmap='jet', drop_large_err=False, rel_err=0.05, **kwargs):
        if ax is None:
            ax = plt.gca()
        if df is None:
            raise ValueError('Empty df')
        df.query('z > 0.0', inplace=True)
        if drop_large_err:
            df.query(f'z_ferr < {rel_err}', inplace=True)

        cmap = copy(plt.cm.get_cmap(cmap))
        cmap.set_under('white')
        kw = dict(cmap=cmap)
        kw.update(kwargs)
        return ax.hist2d(df.x, df.y, weights=df['z'], **kw)

    def rebin1d(self, df, range=(0, 600), bins=30):
        df = df.copy()
        df.query(f'x >= {range[0]} & x <= {range[1]}', inplace=True)

        hist, edges = np.histogram(
            df.x, range=range, bins=bins, weights=df['y'])
        histerr, edges = np.histogram(
            df.x, range=range, bins=bins, weights=df['y_err']**2)
        histerr = np.sqrt(histerr)
        return pd.DataFrame({
            'x': 0.5 * (edges[1:] + edges[:-1]),
            'y': hist,
            'y_err': histerr,
            'y_ferr': np.divide(histerr, hist, where=(hist != 0.0), out=np.zeros_like(histerr))
        })


class reader:
    def __init__(self, path):
        self.path = pathlib.Path(__file__).parent
        self.path = pathlib.Path(self.path, path).resolve()

    def get_names(self, keyword='None'):
        if not self.path.exists():
            raise OSError(f'file not found : {str(self.path)}')

        if self.path.is_file() and os.access(self.path, os.R_OK):
            with root6.TFile(str(self.path)) as file:
                if keyword is None:
                    return [key.GetName() for key in file.GetListOfKeys()]
                else:
                    return [key.GetName() for key in file.GetListOfKeys() if keyword in key.GetName()]

    def get_histogram(self, name=None, keyword=None):
        if name is None:
            names = self.get_names(keyword)
            if len(names) == 0:
                raise ValueError('No objects is returned.')
            name = min(names, key=len)

        with root6.TFile(self.path) as file:
            hist = file.Get(name)
            hist.SetDirectory(0)  # will not be deleted when file is closed
            return hist
