import pandas as pd
import numpy as np
from typing import Literal
import iminuit

from pyamd.e15190 import e15190
from pyamd.utilities import minuit

class Isoscaling:
    """ A class for handling the fitting of the isoscaling relation :math: `R_{21} = C exp(N\\alpha + Z\\beta)`. After setting up the data, i.e. R21 ratio for different particles, the fit is performed by assuming a fixed normalization `C` for all data points. Then, the best normalization is found by iterating over a range of `C` and fitting the data. Finally, the global fit is performed using the best normalization. The result is stored in the attributes `chi2`, `param_names`, `param_values`, `param_errors`. The predicted R21 ratio for a particle with proton Z and neutron N can be calculated by calling `predict` method. This class relies on `iminuit` package for the fitting but also include the equivalent implementation using ROOT.TMinut. The ROOT version is not recommended as it is much slower than the iminuit version.
    """
    def __init__(self):
        self.r21 = dict()
        # minuit attributes to be stored later
        self.is_fitted = False
        self.chi2 = None
        self.param_names = None
        self.param_values = dict()
        self.param_errors = dict()

    def SetR21(self, particle, X):
        self.r21[particle] = X
        
    @staticmethod
    def model(Z, N, norm, alpha, beta):
        return norm * np.exp(N*alpha + Z*beta)

    @staticmethod
    def model_error(Z, N, norm, alpha, beta, norm_err, alpha_err, beta_err):
        return Isoscaling.model(Z, N, norm, alpha, beta) * np.sqrt(N**2 * alpha_err**2 + Z**2 * beta_err**2 + norm_err**2 / norm**2)

    def _fit(self, X, init_norm=1.0, init_alpha=0.5, init_beta=-0.5, **kwargs):
        """ isoscaling fit for one Pt 
        Parameters
        ----------
        X : np.ndarray
            X[0] : [(Z,N)]
            X[1] : y
            X[2] : yerr
        """
        kw = {
            'ncall' : None, 
            'iterate' : 5, 
            'verbose' : 1,
            'limits' : [[0.8,1.2],[None,None],[None,None]],
            'fix' : [True, False, False],
        }
        kw.update(kwargs)

        def cost_fcn(norm, alpha, beta):
            d2 = (X[1] - Isoscaling.model(*(np.array(list(X[0])).transpose()), norm, alpha, beta)) ** 2
            e2 = X[2] ** 2
            return np.sum(np.divide(d2, e2, where=(e2 != 0.), out=np.zeros_like(d2)))
        
        minuit = iminuit.Minuit(cost_fcn, norm=init_norm, alpha=init_alpha, beta=init_beta)

        minuit.limits = kw['limits']
        for i, para in enumerate(minuit.params):
            minuit.fixed[para.name] = kw['fix'][i]

        minuit.print_level = kw['verbose']
        minuit.migrad(kw['ncall'], kw['iterate'])
        minuit.hesse()
        return minuit

    def fit(self, X, method:Literal['iterate','free']='iterate', norm_step_size=0.01, norm_range=(0.7,1.2)):
        """ performa a global isoscaling fit for all data points
        Parameters
        ----------
        X : np.ndarray
            X[0] : [(Z, N)]
            X[1] : [[y]]
            X[2] : [[yerr]]

        method : str
            `iterate` : fix `norm` across all points, repeat the fit for differenet `norm`
            `free` : do not fix `norm`, fit each data independently
        step_size : float
            used when `method=iterate`, step size varying `norm`.
        range : tuple of float
            used when `method=iterate`, range of `norm`.
        """

        if method == 'iterate':
            minimum_chi2 = 1e10
            minimum_norm = norm_range[0]

            for n in np.arange(*norm_range, norm_step_size):
                chi2 = np.sum([
                    self._fit(
                        X = np.array([list(X[0]), X[1][i], X[2][i]], dtype=object),
                        init_norm = n,
                        fix = [True, False, False]
                    ).fval for i in range(X[1].shape[0])
                ])
                if chi2 < minimum_chi2:
                    minimum_chi2 = chi2
                    minimum_norm = n
            
        elif method == 'free':
            minimum_norm = norm_range[0]

        # do the global fit using the best normalization
        minuit_collection = [
            self._fit(
                X = np.array([list(X[0]), X[1][i], X[2][i]], dtype=object),
                init_norm = minimum_norm, 
                fix = [method != 'free', False, False]
            ) for i in range(X[1].shape[0])
        ]

        # saving result
        self.chi2 = np.array(m.fval for m in minuit_collection)
        self.param_names = minuit_collection[0].parameters
        param_values = np.array([m.values for m in minuit_collection]).transpose()
        param_errors = np.array([m.errors for m in minuit_collection]).transpose()
        for i, name in enumerate(self.param_names):
            self.param_values[name] = param_values[i]
            self.param_errors[name] = param_errors[i]
        
        self.fitted = True
        return self

    def predict(self, X):
        """ Calculate isoscaling R21 for a particle with proton Z and neutron N
        Parameters
        ----------
        X : np.ndarray
            X[0] : Z
            X[1] : N
        """
        y = Isoscaling.model(*X, self.param_values['norm'], self.param_values['alpha'], self.param_values['beta'])
        yerr = Isoscaling.model_error(*X, self.param_values['norm'], self.param_values['alpha'], self.param_values['beta'], self.param_errors['norm'], self.param_errors['alpha'], self.param_errors['beta'])
        return y, yerr
    
    def fit_TMinuit(self, particles=None, range=(100, 400), drop_neutron=True, fixed_norm=False, norm=0.7, verbose=0):
        """ Fit Isoscaling using TMinuit
        """
        global fcn
        global subdf

        if particles is None:
            particles = ['p', 'd', 't', '3He', '4He']
        self.range_in_fit = range
        self.particles_in_fit = particles

        data = dict()
        for particle in particles:
            data[particle] = self.r21[particle].copy()
            data[particle].query(
                f'x >= {range[0]} & x <= {range[1]}', inplace=True)
            data[particle][['N', 'Z']] = e15190.Particle(
                particle).N, e15190.Particle(particle).Z

        if drop_neutron and 'n' in data:
            data.pop('n')

        data = pd.concat(list(data.values()), ignore_index=True)

        def fcn(npar, gin, f, par, iflag):
            chisq = 0.
            delta = (
                subdf.y - Isoscaling.model(subdf['Z'], subdf['N'], par[0], par[1], par[2])) / subdf['y_err']
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
            self.fit_TMinuit(particles=particles,
                     range=range, drop_neutron=drop_neutron, fixed_norm=True, norm=x)
            if self.chisq < chisq_:
                norm_ = x
                chisq_ = self.chisq
        return norm_

    def predict_TMinuit(self, particle):
        df = self.r21[particle].copy()
        df.query(
            f'x >= {self.range_in_fit[0]} & x <= {self.range_in_fit[1]}', inplace=True)

        y = Isoscaling.model(e15190.Particle(particle).Z, e15190.Particle(
            particle).N, self.norm, self.alpha, self.beta)
        yerr = Isoscaling.model_error(e15190.Particle(particle).Z, e15190.Particle(
            particle).N, self.norm, self.alpha, self.beta, self.norm_err, self.alpha_err, self.beta_err)
        return pd.DataFrame({
            'x': df.x.to_numpy(),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(
                yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })