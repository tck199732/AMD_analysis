import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from pyamd.utilities import helper, minuit
from pyamd.e15190 import e15190
from pyamd.analysis import r21

helper = helper.helper()


class Isoscaling:
    PARTICLES = ['n', 'p', 'd', 't', '3He', '4He']

    def __init__(self):
        self.R21 = dict()

    def SetR21(self, particle, df):
        if isinstance(df, pd.DataFrame):
            self.R21[particle] = r21.R21(df_R21=df)
        elif isinstance(df, r21.R21):
            self.R21[particle] = df

    def AverageLogR21(self, range=(0., 600.)):
        average_log_r21 = [self.R21[particle].AverageLogR21(
            range=range) for particle in self.PARTICLES]
        avg, err = np.array(average_log_r21).transpose()

        return pd.DataFrame({
            'Z': [e15190.particle(particle).Z for particle in self.PARTICLES],
            'N': [e15190.particle(particle).N for particle in self.PARTICLES],
            'range': [range] * len(self.PARTICLES),
            'avg_log_r21': avg,
            'avg_log_r21_err': err
        }, index=[particle for particle in self.PARTICLES])

    def plotAverageLogR21(self, ax=None, range=(0., 600), groupby='Z', xlabel='N', drop_neutron=False, **kwargs):

        df = self.AverageLogR21(range=range)
        if drop_neutron:
            df.drop('n', axis=0, inplace=True)
        for val, subdf in df.groupby(groupby):
            subdf_ = pd.DataFrame({
                'x': subdf[xlabel],
                'y': subdf['avg_log_r21'],
                'y_err': subdf['avg_log_r21_err']
            })
            helper.plot1d(ax, subdf_, label=f'Z={val}', **kwargs)

            if subdf_.shape[0] == 1:
                continue
            par, cov = curve_fit(
                lambda x, m, c: m*np.array(x) + c,
                subdf_.x.to_numpy(),
                subdf_.y.to_numpy(),
                sigma=subdf_['y_err'].to_numpy(),
                absolute_sigma=True
            )
            ax.errorbar(subdf_.x, par[0] * subdf_.x + par[1], fmt='--')

        return ax

    def model(self, Z, N, norm, alpha, beta):
        return norm * np.exp(N*alpha + Z*beta)

    def model_err(self, Z, N, norm, alpha, beta, norm_err, alpha_err, beta_err):
        return self.model(Z, N, norm, alpha, beta) * np.sqrt(N**2*alpha_err**2 + Z**2*beta_err**2 + norm_err**2/norm**2)

    def fit(self, particles=None, range=(100, 400), drop_neutron=True, fixed_norm=False, norm=0.7, verbose=0):
        global fcn
        global subdf

        if particles is None:
            particles = self.PARTICLES
        self.range_in_fit = range
        self.particles_in_fit = particles

        data = dict()
        for particle in particles:
            data[particle] = self.R21[particle].df.copy()
            data[particle].query(
                f'x >= {range[0]} & x <= {range[1]}', inplace=True)
            data[particle][['N', 'Z']] = e15190.particle(
                particle).N, e15190.particle(particle).Z

        if drop_neutron and 'n' in data:
            data.pop('n')

        data = pd.concat(list(data.values()), ignore_index=True)

        def fcn(npar, gin, f, par, iflag):
            chisq = 0.
            delta = (
                subdf.y - self.model(subdf['Z'], subdf['N'], par[0], par[1], par[2])) / subdf['y_err']
            chisq += np.sum(delta**2)
            f.value = chisq

        self.norm = []
        self.alpha = []
        self.beta = []
        self.norm_err = []
        self.alpha_err = []
        self.beta_err = []

        self.chisq = 0.
        for val, subdf_ in data.groupby('x'):
            subdf = subdf_
            gMinuit = minuit.TMinuit(3, fcn)
            gMinuit.set_parameter(0, 'norm', norm, 0.01, fixed=fixed_norm)
            gMinuit.set_parameter(1, 'alpha', 0.5, 0.01)
            gMinuit.set_parameter(2, 'beta', 0.5, 0.01)
            gMinuit.fit(verbose=verbose)

            # getters = [gMinuit.get_parameter, gMinuit.get_error]
            # names = [gMinuit.get_parameter_name(n) for n in range(3)]
            self.norm.append(gMinuit.get_parameter('norm'))
            self.alpha.append(gMinuit.get_parameter('alpha'))
            self.beta.append(gMinuit.get_parameter('beta'))
            self.norm_err.append(gMinuit.get_error('norm'))
            self.alpha_err.append(gMinuit.get_error('alpha'))
            self.beta_err.append(gMinuit.get_error('beta'))
            self.chisq += gMinuit.get_chisq()

        self.norm = np.array(self.norm)
        self.alpha = np.array(self.alpha)
        self.beta = np.array(self.beta)
        self.norm_err = np.array(self.norm_err)
        self.alpha_err = np.array(self.alpha_err)
        self.beta_err = np.array(self.beta_err)
        return self

    def best_normalization(self, particles=None, range=(100, 400), drop_neutron=True, norm_range=(0.7, 1.2), norm_step=0.01):

        chisq_ = 1e5
        norm_ = norm_range[0]
        for x in np.arange(*norm_range, norm_step):
            self.fit(particles=particles,
                     range=range, drop_neutron=drop_neutron, fixed_norm=True, norm=x)
            if self.chisq < chisq_:
                norm_ = x
                chisq_ = self.chisq

        return norm_

    def predict(self, particle):
        df = self.R21[particle].df.copy()
        df.query(
            f'x >= {self.range_in_fit[0]} & x <= {self.range_in_fit[1]}', inplace=True)

        y = self.model(e15190.particle(particle).Z, e15190.particle(
            particle).N, self.norm, self.alpha, self.beta)
        yerr = self.model_err(e15190.particle(particle).Z, e15190.particle(
            particle).N, self.norm, self.alpha, self.beta, self.norm_err, self.alpha_err, self.beta_err)
        return pd.DataFrame({
            'x': df.x.to_numpy(),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(
                yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })

    def plotIsoscaling(self, ax=None, particle=None, **kwargs):
        if particle is None:
            raise ValueError('particle not given')

        df = self.predict(particle)
        if ax is None:
            ax = plt.gca()
        ax.fill_between(df.x, df.y-df['y_err'], df.y+df['y_err'], **kwargs)
        # ax.errorbar(df.x, df.y, yerr=df['y_err'], **kwargs)
        return ax

    def plotChemPotentialAlpha(self, ax=None, **kwargs):
        df = self.R21['p'].df.copy()
        df.query(
            f'x >= {self.range_in_fit[0]} & x <= {self.range_in_fit[1]}', inplace=True)
        x = df.x.to_numpy()

        df = pd.DataFrame({
            'x': x,
            'y': self.alpha,
            'y_err': self.alpha_err,
            'y_ferr': np.divide(self.alpha_err, self.alpha, where=(self.alpha != 0.0), out=np.zeros_like(self.alpha_err))
        })
        return helper.plot1d(ax, df, **kwargs)

    def plotChemPotentialBeta(self, ax=None, **kwargs):
        df = self.R21['p'].df.copy()
        df.query(
            f'x >= {self.range_in_fit[0]} & x <= {self.range_in_fit[1]}', inplace=True)
        x = df.x.to_numpy()
        df = pd.DataFrame({
            'x': x,
            'y': self.beta,
            'y_err': self.beta_err,
            'y_ferr': np.divide(self.beta_err, self.beta, where=(self.beta != 0.0), out=np.zeros_like(self.beta_err))
        })
        return helper.plot1d(ax, df, drop_zeros=True, **kwargs)
