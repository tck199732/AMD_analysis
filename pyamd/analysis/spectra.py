from pyamd import PROJECT_DIR
from pyamd.e15190 import e15190
from pyamd.utilities import root6, helper
import pathlib
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from copy import copy
import random
from datetime import datetime
random.seed(datetime.now().timestamp())

plot_helper = helper.DataFrameHelper()
hist_reader = root6.HistogramReader()


class SpectraFile:
    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()

    def get_names(self, particle='None'):
        if not self.path.exists():
            raise OSError(f'file not found : {str(self.path)}')

        if self.path.is_file() and os.access(self.path, os.R_OK):
            with root6.TFile(str(self.path)) as file:
                if particle is None:
                    return [key.GetName() for key in file.GetListOfKeys()]
                else:
                    return [key.GetName() for key in file.GetListOfKeys() if key.GetName().endswith('_' + particle)]

    def get_histogram(self, name=None, particle=None, keyword=None):
        if name is None:
            names = self.get_names(particle)
            names = [name for name in names if keyword.lower() in name.lower()]
            if len(names) == 0:
                raise ValueError('No objects is returned.')
            name = min(names, key=len)
            print(f'reading : {name}')

        with root6.TFile(self.path) as file:
            hist = file.Get(name)
            hist.SetDirectory(0)  # will not be deleted when file is closed
            return hist


class PtRapidityLAB:
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

    DIR = f'{str(PROJECT_DIR)}/result/spectra'

    def __init__(self, particle, df=None, path=None, reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(0., 3.), uball_multiplicity=(1, 25), mode='seq', histname=None):

        self.betacms = e15190.reaction(reaction).get_betacms()
        self.beam_rapidity = e15190.reaction(reaction).get_rapidity_beam()

        self.reaction = reaction

        if df is None:
            if path is None:
                path = f'{self.DIR}/{reaction}_{skyrme}_Ncmin{uball_multiplicity[0]}_Ncmax{uball_multiplicity[1]}_b{impact_parameter[1]:.0f}fm.root'

            spectra_file = SpectraFile(path)

            if histname is None:
                # histname in AMD result
                histname = f'h2_pta_rapidity_lab_{mode}_{particle}'

            df = spectra_file.get_histogram(
                particle=particle, keyword=histname)
            df = hist_reader.hist2d_to_df(df, keep_zeros=False)

        if not 'z_ferr' in df.columns:
            df.columns = ['x', 'y', 'z', 'z_err']
            df['z_ferr'] = np.divide(df['z_err'], df['z'], where=(
                df['z'] != 0.0), out=np.zeros_like(df['z_err']))

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

        # In TH2D, projection is done on the bin range (bin1, bin2) instead of numerical range
        # e.g. if bin width is 0.005, then 0.396 is included in ROOT Projection.
        # 0.395 is not included becasue any value exactly the same as the interval would be filled in the bin on the left. so (> and <=)

        df.query(
            f'x > {float(x1)} & x <= {float(x2)} & y > {float(y1)} & y <= {float(y2)}', inplace=True)
        return df

    def correct_coverage(self, df=None, xrange=(0.4, 0.6), yrange=(0., 600.)):
        df = self.query(df, xrange, yrange)
        df.query('z > 0.0', inplace=True)
        df.reset_index(inplace=True, drop=True)
        correction = np.ones(df.shape[0])
        for val, subdf in df.groupby('y'):
            if subdf.shape[0] == 1:
                correction[np.array(subdf.index)] = 0.
            else:
                dx = np.diff(np.array(subdf.x)[0:2])
                x1 = np.min(subdf.x) + dx * np.random.uniform(-0.5, 0.5)
                x2 = np.max(subdf.x) + dx * np.random.uniform(-0.5, 0.5)
                correction[np.array(subdf.index)] = abs(x2-x1)
        df['correction'] = correction

        df.query('correction > 0.0', inplace=True)
        df['correction'] = df['correction'] / np.diff(xrange)
        df['z'] = df['z']/df['correction']
        df['z_err'] = df['z_err']/df['correction']
        df.drop('correction', axis=1, inplace=True)
        return df

    def PtSpectrum(self, xrange=(0.4, 0.6), yrange=(0., 600.), bins=30,  correct_coverage=False, correct_range=(0., 600)):
        if correct_coverage:
            df = self.correct_coverage(xrange=xrange, yrange=correct_range)
        else:
            df = self.query(xrange=xrange, yrange=yrange)
        df.query('z > 0.0', inplace=True)

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
        return plot_helper.plot2d(ax, df, range=range, bins=bins, **kwargs)

    def plotPtSpectrum(self, ax=None, xrange=(0.4, 0.6), yrange=(0, 600), bins=30, correct_coverage=False, correct_range=(0., 600), **kwargs):
        df = self.PtSpectrum(xrange=xrange, yrange=yrange,
                             bins=bins, correct_coverage=correct_coverage, correct_range=correct_range)
        return plot_helper.plot1d(ax, df, **kwargs)

    def plotRapidityCMS(self, ax=None, range=(-0.5, 0.5), bins=100, **kwargs):
        df = self.RapidityCMS(range=range, bins=bins)
        return plot_helper.plot1d(ax, df, **kwargs)

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

    def plotEkinCMS(self, ax=None, range=(0, 200), bins=50, **kwargs):
        df = self.EkinCMS(xrange=range, bins=bins)
        return plot_helper.plot1d(ax, df, **kwargs)

    def plotThetaCMS(self, ax=None, range=(40, 140), bins=100, **kwargs):
        df = self.ThetaCMS(yrange=range, bins=bins)
        return plot_helper.plot1d(ax, df, **kwargs)

    def plotEkinThetaCMS(self, ax=None, range=[[0, 200], [40, 140]], bins=[50, 100], **kwargs):
        df = self.query(xrange=range[0], yrange=range[1])
        return plot_helper.plot2d(ax, df, range=range, bins=bins, **kwargs)


if __name__ == '__main__':

    spectra = PtRapidityLAB(particle='p', reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(
        0., 3.), uball_multiplicity=(11, 25), mode='seq')

    fig, ax = plt.subplots()
    spectra.plotPtRapidity(ax, norm=LogNorm())
    fig.savefig('test.png')
