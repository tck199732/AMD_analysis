import pandas as pd
import numpy as np
from pyamd.utilities import helper
df_helper = helper.DataFrameHelper()
class R21:
    def __init__(self, particle, df_R21 = None, df_prich=None, df_nrich=None):
        self.particle = particle
        if df_R21 is None:
            df_R21 = df_helper.ratio1d(df_prich, df_nrich)
        self.df = df_R21

    def plotR21(self, ax=None, range=(0,600), bins=30, **kwargs):
        df = self.df.copy()
        
        hist, edges = np.histogram(df.x, range=range, bins=bins, weights=df['y'])
        histerr, edges = np.histogram(df.x, range=range, bins=bins, weights=df['y_err']**2)
        histerr = np.sqrt(histerr)
        
        df = pd.DataFrame({
            'x' : 0.5 * (edges[1:]+edges[:-1]),
            'y' : hist, 
            'y_err': histerr,
            'y_ferr' : np.divide(histerr, hist, where=(hist!=0.0), out=np.zeros_like(histerr))
        })
        return df_helper.plot1d(ax, df, drop_zeros=True, drop_large_err=True, rel_err=0.05, **kwargs)
    
    def AveargeR21(self, range=(0., 600.)):
        df = self.df.copy()
        x1 = np.maximum(np.min(df.x), range[0])
        x2 = np.minimum(np.max(df.x), range[1])
        
        df.query(f'y > 0. & y_err > 0.', inplace=True)
        df.query(f'x >= {x1} & x <= {x2}', inplace=True)

        weights = 1./df['y_err']**2
        average = np.average(df.y, weights=weights)
        error = np.sqrt(np.sum((weights * df['y_err'])**2)) / np.sum(weights)

        return average, error

    def AverageLogR21(self,  range=(0., 600.)):
        df = self.df.copy()
        df.query(f'x >= {range[0]} & x <= {range[1]} & y > 0.0 & y_err > 0.0', inplace=True)

        weights = np.where(df['y_ferr'] > 0.0, 1./df['y_ferr']**2, 0.0)
        average = np.average(np.log(df.y), weights=weights)
        error = np.sqrt(
            np.sum((weights * df['y_ferr'])**2)) / np.sum(weights)
        return average, error
    
    def AverageR21Trend(self, range=(100., 400.)):
            
        df = self.df.copy()
        df.query(
            f'y > 0.0 & y_err > 0.0 & x >= {range[0]} & x <= {range[1]}', inplace=True)
        xrange = (np.maximum(np.min(df.x), range[0]), np.minimum(
            np.max(df.x), range[1]))
        x = np.linspace(*xrange, len(df.x))
        y, yerr = np.array([self.AveargeR21(range=(pt, xrange[1])) for pt in x]).transpose()

        return pd.DataFrame({
            'x': x,
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, out=np.zeros_like(yerr), where=(y > 0.0))
        })

    def plotAverageR21Trend(self, ax=None, range=(0., 600.), **kwargs):
        df = self.AverageR21Trend(range=range)
        return df_helper.plot1d(ax, df, **kwargs)

    