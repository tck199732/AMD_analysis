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

# local
from pyamd.utilities import minuit, root6, symwrap, histo

hist_reader = root6.HistogramReader()
# alias
dfhist = histo.histogram

class Multiplicity_ImpactParameter:
    def __init__(self, path, histname='h2_multi_b', **kwargs):
        """ A Class handle calculation from 2D histogram Impact-Parameter VS Multiplicity
        Parameter
        ---------
        path : str, pathlib.Path
            path to the root file containing the histogram of impact-parameter vs multiplicity
        histname : str
            name of the TH2D in the root file
        kwargs : dict
            additional information about the histogram in the root file
        """

        attributes = ['reaction', 'mode', 'table', 'skyrme', 'impact_parameter']
        for attr in attributes:
            setattr(self, attr, kwargs.get(attr, None))

        with root6.TFile(str(path), mode='READ') as file:
            try:
                hist = file[histname]
            except:
                objname = file.keys()[0]
                warnings.warn(
                    f'hist name {histname} not found in root file. We will use the first object in the file, i.e. `{objname}`')
                hist = file[objname]
            
            # number of events used in the histogram (after filtering)
            self.nentries = hist.GetEntries()
            # number of events in simulation (used as normalization factor)
            self.integral = self.nentries / hist.Integral()

            self.bins = (hist.GetNbinsX(), hist.GetNbinsY())
            self.binwidths = (hist.GetXaxis().GetBinWidth(1), hist.GetYaxis().GetBinWidth(1))
            self.xrange = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
            self.yrange = (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())
            self.df = hist_reader.hist2d_to_df(hist)

    def get_multiplicity_distribution(self, range=(-0.5, 79.5), cut=(0,10), bins=None, normalize=True):
        """ Projection to get multiplicity distribution
        Parameters
        ----------
        range : tuple of float
            range of multiplicity to be included in the final spectra
        cut : tuple of float
            cut on impact-parameter
        bins : int
            number of bins, if None, use 1 bin for each multiplicity
        normalize : bool
            if True, normalize by bin width
        Return  
        ------
        pd.DataFrame
            dataframe containing the multiplicity distribution, with columns `x`, `y`, `y_err`, `y_ferr`
        """ 
        bins = int(np.diff(range)[0]) if bins is None else bins
        return dfhist.projection(self.df, cut=cut, range=range, bins=bins, axis='x', normalize=normalize)

    def get_impact_parameter_distribution(self, range=(0., 10.), cut=(-0.5, 79.5), bins=100, normalize=True):
        """ get impact-parameter distribution 
        Parameter
        ---------
        range : tuple of float 
            range of b
        cut : tuple of float
            cut on multiplicity
        bins : int  
            number of bins
        normalize : bool    
            if True, normalize by bin width
        Return
        ------
        pd.DataFrame
            dataframe containing the impact-parameter distribution, with columns `x`, `y`, `y_err`, `y_ferr`
        """
        return dfhist.projection(self.df, cut=cut, range=range, bins=bins, axis='y', normalize=normalize)
    
    def get_dsigma_db(self, range=(0., 10.), cut=(-0.5, 79.5), bins=20, unit=None):
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
        # normalized impact parameter spectra (1 / N_{simulation}) dN/db
        df = self.get_impact_parameter_distribution(
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

        return pd.DataFrame({
            'x': df.x,
            'y': dsigma_db,
            'y_err': error,
            'y_ferr': np.divide(error, dsigma_db, out=np.zeros_like(error), where=(dsigma_db > 0.0))
        })
        
    def get_sigma_b(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=None):
        """ Integrate :math: `d\sigma/db` to get :math: `\sigma_b`
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
        df = self.get_dsigma_db(range, cut, bins, unit)
        db = df['x'][1] - df['x'][0]
        return np.sum(df['y']) * db, np.sqrt(np.sum(df['y_err'])) * db

    @staticmethod
    def _get_b_from_sigma(sigma):
        return np.sqrt(sigma / np.pi)

    def get_bmax(self, range=(0., 10.), cut=(-0.5, 24.5), bins=20, unit=None):
        sig, err = self.get_sigma_b(range, cut, bins, unit)
        return Multiplicity_ImpactParameter._get_b_from_sigma(sig), np.divide(err / 2., np.sqrt(np.pi * sig), out=np.zeros_like(err), where=(sig > 0.0))

    def get_cummulative_multiplicity_distribution(self, ranges=[-0.5,24.5], cut=[0.,10.], bins=None):
        """ Get cummulative multiplicity distribution 
        Parameters
        ----------
        ranges : list of float
            range of multiplicity to be included in the final spectra
        cut : list of float
            cut on impact-parameter                             
        bins : int
            number of bins, if None, use 1 bin for each multiplicity
        Return
        ------
        pd.DataFrame
            dataframe containing the cummulative multiplicity distribution, with columns `x`, `y`, `y_err`, `y_ferr`
        """
        bins = int(np.round(np.diff(ranges))) if bins is None else bins
        df = self.get_multiplicity_distribution(range=ranges, cut=cut, bins=bins)

        
        P = np.cumsum(df['y'][::-1])[::-1]
        P_err = np.sqrt(np.cumsum(df['y_err'][::-1]**2))[::-1]

        P /= P[0]
        P_ferr = np.divide(P_err, P, out=np.zeros_like(P_err), where=(P!=0.))
        P_err  = P * np.sqrt( P_ferr ** 2. + (df['y_err'][0] / df['y'][0])**2)

        return pd.DataFrame({
            'x' : df['x'],
            'y' : P,
            'y_err' : P_err,
            'y_ferr' : np.divide(P_err, P, out=np.zeros_like(P_err), where=(P!=0.))
        })

    def get_bijective_map(self, ranges=[[-0.5, 24.5], [0., 10.]], bins=[None, 20], unit=None):
        """ Calculate the mapping from charged-particle multiplicty to Impact-Parameter
        Parameters
        ----------
        ranges : list of tuple of float
            range of multiplicity and impact-parameter
        bins : list of int  
            number of bins for multiplicity and impact-parameter
        unit : str  
            `None` : use the same unit as in simulation, i.e. fm^2
            `b` : barn  
            `mb` : 0.1 barn
        Return  
        ------
        pd.DataFrame
            dataframe containing the bijective map, with columns `x`, `y`, `y_err`, `y_ferr
        """
        PNC = self.get_cummulative_multiplicity_distribution(ranges=ranges[0], cut=ranges[1], bins=bins[0])

        bmax, bmax_err = self.get_bmax(
            range=ranges[1], cut=ranges[0], bins=bins[1], unit=unit)

        y = bmax * np.sqrt(PNC['y'])
        yerr = np.sqrt(
            PNC['y'] * bmax_err ** 2 + \
            (0.5 * bmax * PNC['y_ferr']) ** 2.
        )

        return pd.DataFrame({
            'x': PNC['x'].to_numpy(),
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, out=np.zeros_like(yerr), where=(y > 0.0))
        }).query('y != 0.0 and y_ferr != 0.0')

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
            df = self.get_dsigma_db(
                range=range, bins=bins, unit=None)

        cost_fcn = LeastSquares(
            df['x'], df['y'], df['y_err'], Multiplicity_ImpactParameter.model_dsigma_db)

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
            'y_err': yerr.value,
            'y_ferr': np.divide(yerr.value, y.value, out=np.zeros_like(yerr.value), where=(y.value > 0.0))
        })

        if output_info:
            return df, fit_info
        return df

    def fit_dsigma_db_TMinuit(self, df=None, range=(0., 10.), bins=20, unit=None):
        """ Fit the differential cross-section using TMinuit. This function is kept for legacy purpose, use `fit_dsigma_db_iminuit` instead.
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
            df = self.get_dsigma_db(
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
            'y_ferr': np.divide(yerr.value, y.value, out=np.zeros_like(yerr.value), where=(y.value > 0.0))
        })