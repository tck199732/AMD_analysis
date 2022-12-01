import pandas as pd
import numpy as np
import iminuit
from pyamd.e15190 import e15190
from pyamd.utilities import helper, minuit
df_helper = helper.DataFrameHelper()


class Isoscaling:
    def __init__(self, reaction1, reaction2):
        self.reaction1 = reaction1
        self.reaction2 = reaction2
        self.r21 = dict()

    def SetR21(self, particle, df):
        self.r21[particle] = df

    def _model(self, Z, N, norm, alpha, beta):
        return norm * np.exp(N*alpha + Z*beta)

    def _model_err(self, Z, N, norm, alpha, beta, norm_err, alpha_err, beta_err):
        return self.model(Z, N, norm, alpha, beta) * np.sqrt(N**2*alpha_err**2 + Z**2*beta_err**2 + norm_err**2/norm**2)

    # def ifit(self, range=(100, 400), drop_neutron=True, fixed_norm=False, norm=0.7, verbose=0):

    #     data = dict()
    #     for i, (pn, df) in enumerate(self.r21.items()):
    #         zid = e15190.particle(pn).Z
    #         nid = e15190.particle(pn).N
    #         data[pn] = df.copy().query(f'x >= {range[0]} & x <= {range[1]}')
    #         data[pn]['N'] = [nid] * data[pn].shape[0]
    #         data[pn]['Z'] = [zid] * data[pn].shape[0]

    #     data = pd.concat(list(data.values()), ignore_index=True)

    #     # models[pn] = lambda pn, norm, alpha, beta: self.model(zid, nid, norm, alpha, beta)

    #     # cost_fcns[pn] = iminuit.cost.LeastSquares(df.x, df.y, df['y_eff'], self.m)

    #     # minuit = iminuit.Minuit(fcn, )

    # fit using TMinuit

    def fit(self, particles=None, range=(100, 400), drop_neutron=True, fixed_norm=False, norm=0.7, verbose=0):
        global fcn
        global subdf

        if particles is None:
            particles = ['p', 'd', 't', '3He', '4He']
        self.range_in_fit = range
        self.particles_in_fit = particles

        data = dict()
        for particle in particles:
            data[particle] = self.r21[particle].df.copy()
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
        df = self.r21[particle].df.copy()
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
        df = self.r21['p'].df.copy()
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
        df = self.r21['p'].df.copy()
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
