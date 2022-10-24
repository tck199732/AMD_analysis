import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pyamd.utilities import root6, helper
from pyamd.e15190 import e15190
from pyamd.analysis import spectra

helper = helper.helper()


class Coalescence:
    def __init__(self, reaction, skyrme=None):
        self.reaction = reaction
        self.skyrme = skyrme
        self.spectra = dict()

    def add(self, particle, df):
        self.spectra[particle] = df

    def pseudo_neutron(self, range=(0, 600), bins=30):
        if not (set(['p', 't', '3He']).issubset(self.spectra)):
            raise ValueError(
                'missing spectra for calculating pseudo neutron.')

        df_p = self.spectra['p'].copy()
        df_t = self.spectra['t'].copy()
        df_3He = self.spectra['3He'].copy()

        helper.rebin1d(df_p, range=range, bins=bins)
        helper.rebin1d(df_t, range=range, bins=bins)
        helper.rebin1d(df_3He, range=range, bins=bins)

        y = df_p.y * df_t.y
        y = np.where((df_3He.y != 0.0), y/df_3He.y, 0.0)
        yerr = y * np.sqrt(df_p['y_ferr']**2 +
                           df_t['y_err']**2 + df_3He['y_err']**2)

        return pd.DataFrame({
            'x': df_p.x.to_numpy(),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })

    def npRatio(self, range=(0, 600), bins=30, pseudo_neutron=False):
        if not 'p' in self.spectra:
            raise ValueError('proton spectrum is absent.')

        if not ('n' in self.spectra) and not pseudo_neutron:
            raise ValueError('No neutron spectra.')

        if pseudo_neutron:
            df1 = self.pseudo_neutron(range=range, bins=bins)
        else:
            df1 = self.spectra['n'].copy()

        df2 = self.spectra['p'].copy()
        df1 = helper.rebin1d(df1, range=range, bins=bins)
        df2 = helper.rebin1d(df2, range=range, bins=bins)
        return helper.ratio1d(df2, df1)

    def plot_npRatio(self, ax=None, range=(0, 600), bins=30, pseudo_neutron=False, **kwargs):
        df = self.npRatio(range=range, bins=bins,
                          pseudo_neutron=pseudo_neutron)
        return helper.plot1d(ax, df, drop_zeros=True, drop_large_err=True, rel_err=0.05, **kwargs)

    def ChemicalTemperature(self, bins=30, range=(0, 600)):

        if not (set(['d', '4He', 't', '3He']).issubset(self.spectra)):
            raise ValueError(
                'missing spectra for calculating chemical temperature.')

        df_d = self.spectra['d'].copy()
        df_4He = self.spectra['4He'].copy()
        df_t = self.spectra['t'].copy()
        df_3He = self.spectra['3He'].copy()

        helper.rebin1d(df_d, range=range, bins=bins)
        helper.rebin1d(df_4He, range=range, bins=bins)
        helper.rebin1d(df_t, range=range, bins=bins)
        helper.rebin1d(df_3He, range=range, bins=bins)

        num = df_d['y'] * df_4He['y']
        den = df_t['y'] * df_3He['y']

        T = np.divide(num, den, where=(
            (num != 0.0) & (den != 0.0)), out=np.zeros_like(num))
        err = T * np.sqrt(df_d['y_ferr']**2 + df_4He['y_ferr']
                          ** 2 + df_t['y_ferr']**2 + df_3He['y_ferr']**2)
        err = err / np.log(T, where=(T > 0.0))**2

        T, err = (
            14.3 / np.log(T, where=(T > 0.0)),
            14.3 * np.divide(err, T, out=np.zeros_like(err), where=(T > 0.0))
        )

        T, bin_edges = np.histogram(
            df_d['x'], bins=bins, range=range, weights=T)
        err, bin_edges = np.histogram(
            df_d['x'], bins=bins, range=range, weights=err**2)
        err = np.sqrt(err)

        x_edges = np.linspace(*range, bins+1)
        return pd.DataFrame({
            'x': 0.5 * (x_edges[1:] + x_edges[:-1]),
            'y': T,
            'y_err': err,
            'y_ferr': np.divide(err, T, out=np.zeros_like(err), where=(T > 0.0))
        })

    def plotChemicalTemperature(self, ax=None, range=(0., 600), bins=30, **kwargs):
        df = self.ChemicalTemperature(bins=bins, range=range)
        return helper.plot1d(ax, df, drop_zeros=True, drop_large_err=True, rel_err=0.05, **kwargs)

    def coalescence(self, particles=None, coal='p', range=(0, 600), bins=30, pseudo_neutron=True):

        if pseudo_neutron:
            if 'n' in self.spectra:
                self.spectra.pop('n')
            self.spectra['n'] = self.pseudo_neutron(range=range, bins=bins)
        if particles is None:
            particles = [particle for particle in self.spectra]

        spectra = [self.spectra[particle] for particle in particles]
        keys = [(e15190.particle(particle).Z, e15190.particle(particle).N)
                for particle in particles]

        spectra = [helper.rebin1d(df, range=range, bins=bins)
                   for df in spectra]

        df = pd.concat(spectra, keys=keys, names=['Z', 'N'])
        if coal == 'p':
            y = df.groupby(['x']).apply(
                lambda x: np.sum(x.y * x.index.get_level_values('Z').astype(float)))
        elif coal == 'n':
            y = df.groupby(['x']).apply(
                lambda x: np.sum(x.y * x.index.get_level_values('N').astype(float)))

        yerr = df.groupby(['x']).apply(
            lambda x: np.sum(x['y_err']**2))
        yerr = np.sqrt(yerr)

        edges = np.linspace(*range, bins+1)
        return pd.DataFrame({
            'x': 0.5 * (edges[1:]+edges[:-1]),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })

    def coalescence_n(self, df_n=None, df_p=None, df_d=None, df_t=None, df_3He=None, df_4He=None, range=(0,600),bins=30):
        df_n = helper.rebin1d(df_n, range=range, bins=bins)
        df_p = helper.rebin1d(df_p, range=range, bins=bins)
        df_d = helper.rebin1d(df_d, range=range, bins=bins)
        df_t = helper.rebin1d(df_t, range=range, bins=bins)
        df_3He = helper.rebin1d(df_3He, range=range, bins=bins)
        df_4He = helper.rebin1d(df_4He, range=range, bins=bins)

        edges = np.linspace(*range, bins+1)
        y =  df_n.y + df_d.y + 2.*df_t.y + df_3He.y + 2.*df_4He.y
        yerr = np.sqrt(df_n['y_err']**2  + df_d['y_err']**2 + 4. * df_t['y_err']**2 + df_3He['y_err']**2 + 4. * df_4He['y_err']**2)
        
        return pd.DataFrame({
            'x': 0.5 * (edges[1:]+edges[:-1]),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        }) 

    def coalescence_p(self, df_n=None, df_p=None, df_d=None, df_t=None, df_3He=None, df_4He=None, range=(0,600),bins=30):
        df_n = helper.rebin1d(df_n, range=range, bins=bins)
        df_p = helper.rebin1d(df_p, range=range, bins=bins)
        df_d = helper.rebin1d(df_d, range=range, bins=bins)
        df_t = helper.rebin1d(df_t, range=range, bins=bins)
        df_3He = helper.rebin1d(df_3He, range=range, bins=bins)
        df_4He = helper.rebin1d(df_4He, range=range, bins=bins)

        edges = np.linspace(*range, bins+1)
        y =  df_p.y + df_d.y + df_t.y + 2. * df_3He.y + 2. * df_4He.y
        yerr = np.sqrt(df_p['y_err']**2  + df_d['y_err']**2 +  df_t['y_err']**2 + 4. * df_3He['y_err']**2 + 4. * df_4He['y_err']**2)
        
        return pd.DataFrame({
            'x': 0.5 * (edges[1:]+edges[:-1]),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        }) 
        


    def plotCoalescence(self, ax=None, df=None, particles=None, coal='p', range=(0, 600), bins=30, pseudo_neutron=True, **kwargs):
        if df is None:
            df = self.coalescence(
                particles, coal=coal, range=range, bins=bins, pseudo_neutron=pseudo_neutron)
        return helper.plot1d(ax, df, **kwargs)

# class DoubleRatio:
#     # npratio(nrich)  / npratio(prich)
#     # coal-npratio(nrich)  / coal-npratio(prich)

#     def __init__(self, particle, df_nrich, df_prich):
#         self.particle = particle
#         self.df = helper.ratio1d(df_prich, df_nrich)

#     def plot_dr(self, ax=None, range=(0,600), bins=30, **kwargs):
#         df = self.df.copy()
#         df = helper.rebin1d(df, range=range, bins=bins)
#         return helper.plot1d(ax, df, **kwargs)


    