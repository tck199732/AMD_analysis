import numpy as np
import pandas as pd

class histogram_handler:
    """
    A class for handling 1D / 2D histograms in terms of pd.DataFrame. By default, a dataframe representing a 2D histogram should have the following columns:
    - x : the x-axis
    - y : the y-axis
    - y_err : the error of y
    - y_ferr : the fractional error of y
    - z : the z-axis
    - z_err : the error of z
    - z_ferr : the fractional error of z
    fractions of the error of y and z are optional.
    """
    @staticmethod
    def projection(df, cut=None, range=None, bins=None, axis='x', normalize=True):
        """ Project a 2D histogram (in terms of pd.DataFrame) along an axis
        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the histogram to be projected. 
        cut : tuple, optional
            Range of the axis to be cut. If None, no cut is applied.
        range : tuple, optional
            Range of the axis to be projected. If None, the range of the dataframe is used.
        bins : int, optional
            Number of bins to be used. If None, the number of bins of the dataframe is used.
        axis : str, optional
            Axis along which the projection is performed. Must be 'x' or 'y'.
        normalize : bool, optional
            If True, the projection is normalized to the bin width.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the projected histogram.
        """

        df = df.copy()
        if df.shape[1] < 4:
            raise ValueError('Dataframe must have at least 3 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'z', 'z_err']
        elif df.shape[1] == 5:
            df.columns = ['x', 'y', 'z', 'z_err', 'z_ferr']
        else:
            raise ValueError('Dataframe must have at most 5 columns')
        
        cut_axis = 'y' if axis == 'x' else 'x'
        if cut is not None:
            bw = df[cut_axis][1] - df[cut_axis][0]
            cut = (cut[0] - bw/2, cut[1] + bw/2)
            df.query(f'{cut_axis} > {cut[0]} & {cut_axis} < {cut[1]}', inplace=True)

        range = (df[axis].min(), df[axis].max()) if range is None else range
        bins = df.shape[0] if bins is None else bins

        hist, edges = np.histogram(
            df[axis], 
            weights=df['z'],
            range=range, 
            bins=bins
        )

        histerr, _ = np.histogram(
            df[axis],
            weights=df['z_err']**2, 
            range=range,
            bins=bins
        )
        histerr = np.sqrt(histerr)
        norm = (np.diff(range) / bins)[0] if normalize else 1.
        
        return pd.DataFrame({
            'x' : (edges[1:] + edges[:-1])/2,
            'y' : hist / norm,
            'y_err' : histerr / norm,
            'y_ferr' : np.divide(histerr, hist, out=np.zeros_like(hist), where=hist!=0)
        })

    @staticmethod
    def profile(df, cut=None, axis='x', option=None, scale=1.):
        """
        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the 2D histogram to be profiled.
        cut : tuple, optional
            Range of the axis to be cut. If None, no cut is applied.
        axis : str, optional
            Axis along which the profile is performed. Must be 'x' or 'y'.
        option : str, optional
            option for the error calculation in the same way as in ROOT, i.e. 
            `None` : std. dev. divided by sqrt(scale), where scale is the number of entries in the bin 
            `'s'` : std. dev.
            `'g'` : 1./sqrt(W) where W is the sum of the weights
        scale : float, optional
            Scale factor for the error if option is None.
        """
        df = df.copy()
        if df.shape[1] < 4:
            raise ValueError('Dataframe must have at least 3 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'z', 'z_err']
        elif df.shape[1] == 5:
            df.columns = ['x', 'y', 'z', 'z_err', 'z_ferr']
        else:
            raise ValueError('Dataframe must have at most 5 columns')
        
        df.query('z > 0 & z_err > 0', inplace=True)

        cut_axis = 'y' if axis == 'x' else 'x'
        if cut is not None:
            bw = df[cut_axis][1] - df[cut_axis][0]
            cut = (cut[0] - bw/2, cut[1] + bw/2)
            df.query(f'{cut_axis} > {cut[0]} & {cut_axis} < {cut[1]}', inplace=True)

        hist = df.groupby(cut_axis).apply(
            lambda subdf: np.average(subdf[axis], weights=subdf['z'])
        )

        x = hist.index.to_numpy()
        hist = hist.to_numpy()

        if option is None:

            entries = df.groupby(cut_axis).apply(
                lambda subdf: np.sum(subdf['z'])
            )
            entries *= scale
            mean2 = df.groupby(cut_axis).apply(
                lambda subdf: np.average(subdf[axis] ** 2, weights=subdf['z'])
            )
            histerr = np.sqrt(mean2 - hist ** 2) / np.sqrt(entries)

        elif option == 's':
            mean2 = df.groupby(cut_axis).apply(
                lambda subdf: np.average(subdf[axis] ** 2, weights=subdf['z'])
            )
            histerr = np.sqrt(mean2 - hist ** 2)

        elif option == 'g':
            histerr = df.groupby(cut_axis).apply(
                lambda subdf: 1. / np.sqrt(np.sum(subdf['z_err']))
            )
        
        histerr = histerr.to_numpy()
        return pd.DataFrame({
            'x' : x,
            'y' : hist,
            'y_err' : histerr,
            'y_ferr' : np.divide(histerr, hist, out=np.zeros_like(hist), where=hist!=0)
        })

    @staticmethod
    def query(df, xrange, yrange=None):
        df = df.copy()
        xbw = np.diff(xrange)[0] / df.shape[0]
        xrange = (xrange[0] - xbw/2, xrange[1] + xbw/2)
        df.query(f'x > {xrange[0]} & x < {xrange[1]}', inplace=True)
        if yrange is not None:
            ybw = np.diff(yrange)[0] / df.shape[1]
            yrange = (yrange[0] - ybw/2, yrange[1] + ybw/2)
            df.query(f'y > {yrange[0]} & y < {yrange[1]}', inplace=True)
        return df

    @staticmethod
    def rebin1D(df, range, bins=None, normalize=True):
        df = df.copy()
        if df.shape[1] > 4:
            raise ValueError('Dataframe must have at most 4 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'y_err', 'y_ferr']
        elif df.shape[1] == 3:
            df.columns = ['x', 'y', 'y_err']
        elif df.shape[1] == 2:
            df.columns = ['x', 'y']
            df['y_err'] = 0
        else:
            raise ValueError('Dataframe must have at least 2 columns')
        
        bins = df.shape[0] if bins is None else bins
        hist, edges = np.histogram(
            df['x'],
            weights=df['y'],
            range=range,
            bins=bins
        )

        histerr, _ = np.histogram(
            df['x'],
            weights=df['y_err']**2,
            range=range,
            bins=bins
        )
        histerr = np.sqrt(histerr)
        norm = (np.diff(range) / bins)[0] if normalize else 1.

        return pd.DataFrame({
            'x' : (edges[1:] + edges[:-1])/2,
            'y' : hist / norm,
            'y_err' : histerr / norm,
            'y_ferr' : np.divide(histerr, hist, out=np.zeros_like(hist), where=hist!=0)
        })

    @staticmethod
    def rebin2D(df, range, bins=None, normalize=True):
        df = df.copy()
        if df.shape[1] < 4:
            raise ValueError('Dataframe must have at least 3 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'z', 'z_err']
        elif df.shape[1] == 5:
            df.columns = ['x', 'y', 'z', 'z_err', 'z_ferr']
        else:
            raise ValueError('Dataframe must have at most 5 columns')
        
        if np.array(range).shape != (2, 2):
            raise ValueError('Range must be a 2x2 array')
        if bins is None:
            bins = [df.shape[0], df.shape[1]]
        
        x = np.unique(df['x'])
        y = np.unique(df['y'])
        xbw = x[1] - x[0]
        ybw = y[1] - y[0]
        
        hist, xedges, yedges = np.histogram2d(
            df['x'],
            df['y'],
            weights=df['z'],
            range=range,
            bins=bins
        )

        histerr, _, _ = np.histogram2d(
            df['x'],
            df['y'],
            weights=df['z_err']**2,
            range=range,
            bins=bins
        )
        
        x = (xedges[1:] + xedges[:-1]) / 2.
        y = (yedges[1:] + yedges[:-1]) / 2.

        norm = (x[1] - x[0]) / xbw * (y[1] - y[0]) / ybw if normalize else 1.
        hist /= norm
        histerr = np.sqrt(histerr) / norm

        xx, yy = np.meshgrid(x, y, indexing='ij')
        return pd.DataFrame({
            'x' : xx.flatten(),
            'y' : yy.flatten(),
            'z' : hist.flatten(),
            'z_err' : histerr.flatten(),
            'z_ferr' : np.divide(histerr.flatten(), hist.flatten(), out=np.zeros_like(hist.flatten()), where=hist.flatten()!=0)
        })

    @staticmethod
    def mul(a, b, axis='y'):
        """ Multiply two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely multiply multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional    
            Axis along which the multiplication is performed, depending on the name of the columns. Default is 'y'.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the product of the two histograms.
        """
        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')

        return pd.DataFrame({
            'x' : a['x'],
            axis : a[axis] * b[axis],
            f'{axis}_err' : np.sqrt((a[f'{axis}_err'] * b[axis])**2 + (a[axis] * b[f'{axis}_err'])**2),
            f'{axis}_ferr' : np.sqrt((a[f'{axis}_ferr'] * b[axis])**2 + (a[axis] * b[f'{axis}_ferr'])**2)
        })

    @staticmethod
    def div(a, b, axis='y'):
        """ Divide two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely divide multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional
            Axis along which the division is performed, depending on the name of the columns. Default is 'y'.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the quotient of the two histograms.
        """
        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')
        
        if not 'y_ferr' in a.columns:
            a['y_ferr'] = np.divide(a['y_err'], a['y'], out=np.zeros_like(a['y']), where=a['y']!=0)
        if not 'y_ferr' in b.columns:
            b['y_ferr'] = np.divide(b['y_err'], b['y'], out=np.zeros_like(b['y']), where=b['y']!=0)

        y = np.divide(a[axis], b[axis], out=np.zeros_like(a[axis]), where=b[axis]!=0)
        y_err = y * np.sqrt(a[f'{axis}_ferr'] ** 2. + b[f'{axis}_ferr'] ** 2.)
        return pd.DataFrame({
            'x' : a['x'],
            axis : y,
            f'{axis}_err' : y_err,
            f'{axis}_ferr' : np.divide(y_err, y, out=np.zeros_like(y), where=y!=0)
        })

    @staticmethod
    def sum(a, b, axis='y', sign=1.):
        """ Get the sum of two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely add multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional
            Axis along which the addition is performed, depending on the name of the columns. Default is 'y'.
        sign : float, optional
            Sign of the second histogram. Default is 1. 
        Returns
        -------
        pd.DataFrame
            Dataframe containing the sum of the two histograms.
        """

        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')
        
        y = a[axis] + b[axis] * sign
        y_err = np.sqrt(a[f'{axis}_err']**2 + b[f'{axis}_err']**2)
        return pd.DataFrame({
            'x' : a['x'],
            axis : y,
            f'{axis}_err' : y_err,
            f'{axis}_ferr' : np.divide(y_err, y, out=np.zeros_like(y), where=y!=0)
        })
    
    @staticmethod
    def add(a, b, axis='y'):
        return histogram_handler.sum(a, b, axis, sign=1.)
    
    @staticmethod
    def sub(a, b, axis='y'):
        return histogram_handler.sum(a, b, axis, sign=-1.)
