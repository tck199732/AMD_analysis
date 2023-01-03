import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import copy


class DataFrameHelper:
    def __init__(self):
        pass

    def check_error_information(self, df, label_err='y_err', exception=True):
        has_error = label_err in df.columns
        if not has_error:
            if exception:
                raise ValueError('input does not contain error information.')
            else:
                return False
        return True

    def check_same_axis(self, df1, df2, axis='x', exception=True):
        x1 = df1[axis].to_numpy()
        x2 = df2[axis].to_numpy()
        same_xaxis = np.array_equal(x1, x2)
        if not same_xaxis:
            if exception:
                raise ValueError('numerator and denominator have differnet x.')
            else:
                return False
        return True

    def calculate_relative_error(self, df, label_err='y_err', label_data='y', label='y_ferr', inplace=True):
        if not inplace:
            df = df.copy()
        self.check_error_information(df)
        yferr = np.divide(df[label_err], df[label_data], out=np.zeros_like(
            df[label_err]), where=(df[label_err] != 0))
        df[label] = yferr
        return yferr

    def ratio1d(self, num, den, calculate_error=True):
        """
        Parameters
        ----------
        num : pd.DataFrame
        den : pd.DataFrame
        """
        df1 = num.copy()
        df2 = den.copy()
        self.check_same_axis(df1, df2)

        x = df1.x.to_numpy()
        y = np.where((df1.y != 0.0), df2.y / df1.y, 0.0)

        if calculate_error:
            if not ('y_err' in df1.columns and 'y_err' in df2.columns):
                raise ValueError('input does not contain error information.')

            self.check_error_information(df1)
            self.check_error_information(df2)

            if not 'y_ferr' in df1.columns:
                self.calculate_relative_error(df1)

            if not 'y_ferr' in df2.columns:
                self.calculate_relative_error(df2)

            yerr = y * np.sqrt(df1['y_ferr']**2 + df2['y_ferr']**2)

        return pd.DataFrame({
            'x': x,
            'y': y,
            'y_err': yerr,
            'y_ferr': np.divide(yerr, y, where=(y != 0.0), out=np.zeros_like(yerr))
        })

    def rebin1d(self, df, range=(0, 600), bins=30, normalize=True, inplace=True):

        if not inplace:
            # make df local 
            df = df.copy()

        df.query(f'x >= {range[0]} & x <= {range[1]}', inplace=True)

        hist, edges = np.histogram(
            df.x, range=range, bins=bins, weights=df['y'])
        histerr, edges = np.histogram(
            df.x, range=range, bins=bins, weights=df['y_err']**2)
        histerr = np.sqrt(histerr)

        df.drop(df.index, inplace=True)
        df['x'] = 0.5 * (edges[1:] + edges[:-1])
        df['y'] = hist
        df['y_err'] = histerr
        self.calculate_relative_error(df, inplace=True)

        if normalize:
            df['y'] /= (np.diff(range) / bins)[0]
            df['y_err'] /= (np.diff(range) / bins)[0]

        return df

    def add1d(self, df1, df2, calculate_error=True):
        df1 = df1.copy()
        df2 = df2.copy()
        self.check_same_axis(df1, df2)

        df = pd.DataFrame({
            'x': df1.x,
            'y': df1.y + df2.y
        })

        if calculate_error:
            self.check_error_information(df1)
            self.check_error_information(df2)

            if not 'y_ferr' in df1.columns:
                self.calculate_relative_error(df1)

            if not 'y_ferr' in df2.columns:
                self.calculate_relative_error(df2)

            yerr = np.sqrt(df1['y_err']**2 + df2['y_err']**2)
            df['y_err'] = yerr
            self.calculate_relative_error(df)
        return df

    def average1d(self, df, range):
        df = df.copy()
        weights = np.where(
            (df.x >= range[0]) & (df.x <= range[1]),
            1.,
            0.
        )
        return np.average(df.y, weights=weights)

    def df2text(self, df, filename, columns=['x', 'y', 'y_err', 'y_ferr']):
        df = df.copy()
        df.columns = ['x', 'y', 'y_err', 'y_ferr']
        df.query('y != 0.0', inplace=True)
        df.columns = columns
        df.to_csv(filename, sep=' ', header=True,
                  columns=columns[:-1], index=False, float_format='%.3f')

    def plot1d(self, ax=None, df=None, drop_zeros=True, drop_large_err=False, rel_err=0.05, label_data='y', label_err='y_err', label_ferr='y_ferr', **kwargs):

        if ax is None:
            ax = plt.gca()
        if df is None:
            raise ValueError('Empty dataframe.')

        dfcopy = df.copy()
        if drop_zeros:
            dfcopy.query(f'{label_data} != 0.0', inplace=True)
        if drop_large_err:
            dfcopy.query(f'{label_ferr} < {rel_err}', inplace=True)

        kw = dict(
            fmt='.'
        )
        kw.update(kwargs)
        return ax.errorbar(dfcopy.x, dfcopy.y, yerr=dfcopy[label_err], **kw)

    def plot2d(self, ax=None, df=None, cmap='jet', drop_large_err=False, rel_err=0.05, zlabel='z', label_err='z_err', label_ferr='z_ferr', **kwargs):
        if ax is None:
            ax = plt.gca()
        if df is None:
            raise ValueError('Empty df')

        dfcopy = df.copy()
        dfcopy.query(f'{zlabel} > 0.0', inplace=True)
        if drop_large_err:
            self.check_error_information(dfcopy, label_err=label_err)
            if not label_ferr in dfcopy.columns:
                self.calculate_relative_error(
                    dfcopy, label_err=label_err, label_ferr='z_ferr', label_data='z')
            dfcopy.query(f'{label_ferr} < {rel_err}', inplace=True)

        cmap = copy(plt.cm.get_cmap(cmap))
        cmap.set_under('white')
        kw = dict(cmap=cmap)
        kw.update(kwargs)
        return ax.hist2d(df.x, df.y, weights=df[zlabel], **kw)

    def rebin3d(self, df, bins, ranges=None, inplace=True):

        if not inplace:
            # make df local
            df = df.copy()

        columns = ['x', 'y', 'z', 'content', 'error']
        df.columns = columns
        if ranges is None:
            ranges = [
                (df['x'].min(), df['x'].max()),
                (df['y'].min(), df['y'].max()),
                (df['z'].min(), df['z'].max())
            ]
        edges = [np.linspace(*(ranges[i]), bins[i]+1)
                 for i in range(len(bins))]
        data = df[['x', 'y', 'z']].values
        weights = df['content']
        errors = df['error']
        hist, edges = np.histogramdd(data, bins=edges, weights=weights)
        histerr, edges = np.histogramdd(data, bins=edges, weights=errors**2.)
        histerr = np.sqrt(histerr)

        x, y, z = [
            0.5 * (edge[:-1] + edge[1:]) for edge in edges]
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        
        df.drop(df.index, inplace=True)
        df['x'] = xx.flatten()
        df['y'] = yy.flatten()
        df['z'] = zz.flatten()
        df['content'] = hist.flatten()
        df['error'] = histerr.flatten()

        return df

    def ratio3d(self, num, den, calculate_error=True):

        # check axis here

        num = num.copy()
        den = den.copy()
        columns = ['x', 'y', 'z', 'content', 'error']
        num.columns = columns
        den.columns = columns

        content = np.where((den['content'] != 0.0),
                           num['content'] / den['content'], 0.0)

        if calculate_error:
            # num['ferror'] = np.divide(
            #     num['error'], num['content'], where=(num['content'] > 0.))
            # den['ferror'] = np.divide(
            #     den['error'], den['content'], where=(den['content'] > 0.))
            # error = content * np.sqrt(num['ferror']**2 + den['ferror']**2)
            num['ferror'] = np.where(num['content'] != 0., np.abs(num['error'] / num['content']), 0.)
            den['ferror'] = np.where(den['content'] != 0., np.abs(den['error'] / den['content']), 0.)
            error = content * np.sqrt(num['ferror']**2 + den['ferror']**2)

        return pd.DataFrame({
            'x': num.x,
            'y': num.y,
            'z': num.z,
            'content': content,
            'error': error
        })

    def projection1d(self, df, axis='x', wname=None, ename=None, bins=None, range=None, normalize=False):

        """Projection to 1D histogram represented by pd.DataFrame
        Parameter
        ---------
        df : pd.DataFrame
            N-Dimensional histogram represented by pd.DataFrame
        axis : str
            axis to project, default `x`
        wname : str
            column name of the weight. For 2D histogram, `wname = z`. For 3D histogram, `wname = content`.
        ename : str
            column name of the error.  For 2D histogram, `wname = z_err`. For 3D histogram, `wname = error`.
        normalize : bool
            Projection to 1 axis means summing over the other two axes which is correct for counting histogram. If the histogram represents some function, for instance CF, we should take the average.
        """

        if not axis in df.columns:
            raise ValueError(f'dataframe does not contain column {axis}.')

        if wname is None:
            wname = 'z' if df.shape[1] <= 4 else 'content'
        if ename is None:
            ename = 'z_err' if df.shape[1] <= 4 else 'error'

        if bins is None:
            bins = len(np.unique(df[axis]))

        if range is None:
            range = (np.min(df[axis]), np.max(df[axis]))


        hist, edges = np.histogram(df[axis], bins=bins, range=range, weights=df[wname])

        if ename in df.columns:
            histerr, _ = np.histogram(df[axis], bins=bins, range=range, weights=df[ename]**2)
            histerr = np.sqrt(histerr)
        else:
            histerr = np.zeros(shape=hist.shape)
        
        if normalize:
            norm = np.array([subdf.shape[0] for _, subdf in df.groupby(axis, group_keys=False)])
            hist /= norm
            histerr /= norm


        result = pd.DataFrame({
            'x' : 0.5 * (edges[1:] + edges[:-1]),
            'y' : hist, 
            'y_err' : histerr
        })        
        self.calculate_relative_error(result)

        return result
        
        