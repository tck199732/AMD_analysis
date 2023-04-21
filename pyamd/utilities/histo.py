import numpy as np
import pandas as pd

class histogram:
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
            df.query(f'{cut_axis} >= {cut[0]} & {cut_axis} < {cut[1]}', inplace=True)

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

