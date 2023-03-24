import os
import pathlib
import sys
import subprocess
import numpy as np
import pandas as pd
from typing import Literal
import uproot
import string
import random


try:
    import ROOT
except:
    path_root = subprocess.check_output('which root', shell=True)
    path_root = str(path_root, encoding='utf-8').strip()
    path_root = pathlib.Path(str(path_root)).parent.parent
    path_root = pathlib.Path(path_root, 'lib')
    sys.path.append(str(path_root))
    import ROOT


class RandomNameGenerator:
    def __init__(self):
        self.used_name = set()

    def generate(self, pattern=None, len=10):
        """generate a random string of length 10.
        Parameter
        ---------
        pattern:
            default is string starts with character followed by 9 char/digit. some
             - 'char' : character only
             - 'digit' : digits only
        len : int
            length of the sequence
        """
        while True:
            key = self._generate(pattern=pattern, len=len)
            if not key in self.used_name:
                self.used_name.add(key)
                break
        return key

    def _generate(self, pattern=None, len=10):

        if pattern == 'char':
            entries = [random.choice(string.ascii_letters) for _ in range(len)]
        elif pattern == 'digit':
            entries = [random.choice(string.digits) for _ in range(len)]
        else:
            entries = []
            entries.append(random.choice(string.ascii_letters))
            for _ in range(1, len):
                entries.append(random.choice(
                    string.ascii_letters + string.digits))

        key = ''.join(entries)
        return key


rdnkey_generator = RandomNameGenerator()

class TFile:
    def __init__(self, path, mode='READ', lib : Literal['ROOT', 'uproot'] = 'ROOT'):
        """Constructer of TFile
        Parameters
        ----------
        path : str or pathlib.Path
        mode : str
            READ, RECREATE, etc,... only used when lib == "ROOT".
        lib : str
            - ROOT : use original ROOT TFile
            - uproot : use uproot to read file
        Example
        -------
        with TFile(str(path)) as f:
            h1 = f[key]
        """
        self.path = pathlib.Path(path)
        self.mode = mode.upper()
        self.lib = lib

    def __enter__(self):
        if self.lib == 'ROOT':
            self.file = ROOT.TFile(str(self.path), self.mode)
        elif self.lib == 'uproot':
            self.file = uproot.open(str(self.path))
        return self

    def __getitem__(self, key):
        rdnkey = rdnkey_generator.generate()
        if self.lib == 'ROOT':
            obj = self.file.Get(key)
            obj.SetName(rdnkey)
        if self.lib == 'uproot':
            obj = self.file[key]
            obj.fName = key
        return obj

    def __exit__(self, *args):
        if self.lib == 'ROOT':
            self.file.Close()

    def keys(self):
        if self.lib == 'uproot':
            return self.file.keys()
        if self.lib == 'ROOT':
            return [key.GetName() for key in self.file.GetListOfKeys()]


class HistogramReader:
    def __init__(self):
        pass

    def _ROOT_TH1D_to_df(self, hist,  xname='x', yname='y', keep_zeros=True):
        df = dict()
        columns = [xname, yname, f'{yname}_err']
        getters = [hist.GetXaxis().GetBinCenter, hist.GetBinContent,
                   hist.GetBinError]
        for col, get in zip(columns, getters):
            df[col] = [get(b) for b in range(1, hist.GetNbinsX() + 1)]

        df = pd.DataFrame(df)
        condition = (df[yname] != 0.0)
        df[f'{yname}_ferr'] = np.where(
            condition, np.abs(df[f'{yname}_err']/df[yname]), 0.0)

        return df if keep_zeros else df.query(f'{yname} != 0.0').reset_index(drop=True)

    def _uproot_TH1D_to_df(self, hist,  xname='x', yname='y', keep_zeros=True):
        x = hist.axis().edges()
        x = 0.5 * (x[:-1] + x[1:])
        y = hist.values()
        yerr = hist.errors()
        return pd.DataFrame({
            xname: x,
            yname: y,
            f'{yname}_err': yerr,
            f'{yname}_ferr': np.divide(yerr, y, where=(y != 0.))
        })

    def hist1d_to_df(self, hist, xname='x', yname='y', keep_zeros=True):
        """Convert ROOT.TH1D to pd.DataFrame. Conversion methods depends on using ROOT or uproot.
        Parameter
        ---------
        hist : ROOT.TH1D or uproot.models.TH.Model_TH1D_v{3,4}
            input 1D histogram
        xname, yname : str
            name of the columns
        keep_zeros : bool
            if False, remove all entries with `y == 0`.
        """
        if isinstance(hist, ROOT.TH1D):
            return self._ROOT_TH1D_to_df(hist, xname=xname, yname=yname, keep_zeros=keep_zeros)

        # elif isinstance(hist, uproot.models.TH.Model_TH1D_v3):
        # in case an update of uproot changes the TH1D version, check the base class.
        elif isinstance(hist, uproot.behaviors.TH1.TH1):
            return self._uproot_TH1D_to_df(hist, xname=xname, yname=yname, keep_zeros=keep_zeros)

    def _ROOT_TH2D_to_df(self, histo, xname='x', yname='y', zname='z', keep_zeros=True):
        x = np.array([histo.GetXaxis().GetBinCenter(b)
                     for b in range(1, histo.GetNbinsX() + 1)])
        y = np.array([histo.GetYaxis().GetBinCenter(b)
                     for b in range(1, histo.GetNbinsY() + 1)])

        content = np.array(histo)
        content = content.reshape(len(x) + 2, len(y) + 2, order='F')
        content = content[1:-1, 1:-1]

        error = np.array([histo.GetBinError(b)
                         for b in range((len(x) + 2) * (len(y) + 2))])
        error = error.reshape(len(x) + 2, len(y) + 2, order='F')
        error = error[1:-1, 1:-1]

        xx, yy = np.meshgrid(x, y, indexing='ij')
        df = pd.DataFrame({
            xname: xx.flatten(),
            yname: yy.flatten(),
            zname: content.flatten(),
            f'{zname}_err': error.flatten(),
        })
        mask = (df[zname] != 0.0)
        df[f'{zname}_ferr'] = np.where(
            mask,
            np.abs(df[f'{zname}_err'] / df[zname]),
            0.0
        )
        return df if keep_zeros else df.query(f'{zname} != 0.0').reset_index(drop=True)

    def _uproot_TH2D_to_df(self, histo, xname='x', yname='y', zname='z', keep_zeros=True):
        x = histo.axis('x').edges()
        y = histo.axis('y').edges()
        x = 0.5 * (x[:-1] + x[1:])
        y = 0.5 * (y[:-1] + y[1:])
        xx, yy = np.meshgrid(x, y, indexing='ij')
        z = histo.values().flatten()
        zerr = histo.errors().flatten()

        df = pd.DataFrame({
            xname: xx.flatten(),
            yname: yy.flatten(),
            zname: z,
            f'{zname}_err': zerr
        })
        mask = (df[zname] != 0.0)
        df[f'{zname}_ferr'] = np.where(
            mask,
            np.abs(df[f'{zname}_err'] / df[zname]),
            0.0
        )
        return df if keep_zeros else df.query(f'{zname} != 0.0').reset_index(drop=True)

    def hist2d_to_df(self, histo, xname='x', yname='y', zname='z', keep_zeros=True):
        """Convert ROOT.TH2D to pd.DataFrame.
        Parameter
        ---------
        histo : ROOT.TH2D or uproot.behaviors.TH2.TH2
        xname, yname, zname : str
            name of the pd.DataFrame columns
        Return
        ------
        pandas.DataFrame

        Example
        -------
        # for ROOT.TH2D
        h2 = ROOT.TH2D('test', '', 100, -1., 1., 100, -1., 1.)
        for _ in range(100):
            h2.Fill(ROOT.gRandom.Gaus(0.), ROOT.gRandom.Gaus(0.5))
        df = hist_reader.hist2d_to_df(h2)

        # for uproot.behaviors.TH2.TH2
        h2 = uproot.open(str(path))[hname]
        df = hist_reader.hist2d_to_df(h2)
        """
        if isinstance(histo, ROOT.TH2D):
            return self._ROOT_TH2D_to_df(histo, xname, yname, zname, keep_zeros)

        if isinstance(histo, uproot.behaviors.TH2.TH2):
            return self._uproot_TH2D_to_df(histo, xname, yname, zname, keep_zeros)

    def _uproot_TH3D_to_df(self, histo, xname='x', yname='y', zname='z', wname='content', ename='error', keep_zeros=False):
        xedges = histo.axis('x').edges()
        yedges = histo.axis('y').edges()
        zedges = histo.axis('z').edges()

        x = 0.5 * (xedges[1:] + xedges[:-1])
        y = 0.5 * (yedges[1:] + yedges[:-1])
        z = 0.5 * (zedges[1:] + zedges[:-1])

        content = histo.values()
        errors = histo.errors()
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

        df = pd.DataFrame({
            xname: xx.flatten(),
            yname: yy.flatten(),
            zname: zz.flatten(),
            wname: content.flatten(),
            ename: errors.flatten()
        })

        return df.query('content > 0.') if keep_zeros else df

    def _ROOT_TH3D_to_df(self, histo, xname='x', yname='y', zname='z', wname='content', ename='error', keep_zeros=False):
        x = np.array([histo.GetXaxis().GetBinCenter(b)
                     for b in range(1, histo.GetNbinsX() + 1)])
        y = np.array([histo.GetYaxis().GetBinCenter(b)
                     for b in range(1, histo.GetNbinsY() + 1)])
        z = np.array([histo.GetZaxis().GetBinCenter(b)
                     for b in range(1, histo.GetNbinsZ() + 1)])
        
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

        content = np.array(histo).reshape(len(x)+2, len(y)+2, len(z)+2, order='F')

        content = content[1:-1, 1:-1, 1:-1]

        error = np.array([histo.GetBinError(b)
                        for b in range((len(x) + 2) * (len(y) + 2) * (len(z)+2))])
        error = error.reshape(len(x) + 2, len(y) + 2, len(z)+2, order='F')
        error = error[1:-1, 1:-1, 1:-1]

        df = pd.DataFrame({
            xname: xx.flatten(),
            yname: yy.flatten(),
            zname: zz.flatten(),
            wname: content.flatten(),
            ename: error.flatten()
        })
        mask = (df[wname] != 0.0)
        df[f'f{ename}'] = np.where(
            mask,
            np.abs(df[ename] / df[wname]),
            0.0
        )
        return df if keep_zeros else df.query(f'{wname} != 0.0').reset_index(drop=True)

    def hist3d_to_df(self, histo, xname='x', yname='y', zname='z', wname='content', ename='error', keep_zeros=False):
        """convert a uproot.TH3D to pd.DataFrame
        histo : uproot.models.TH.Model_TH3D_v4
            obtained from `uproot.open(path)[hname]`
        """
        if isinstance(histo, ROOT.TH3D):
            return self._ROOT_TH3D_to_df(histo, xname, yname, zname, wname, ename, keep_zeros)

        if isinstance(histo, uproot.behaviors.TH3.TH3):
            return self._uproot_TH3D_to_df(histo, xname, yname, zname, wname, ename, keep_zeros)

    def projection(self, hist, name=None, axis:Literal['x','y']='x', range=(0, 1e10)):
        """ Project a 2D histogram to 1D
        Parameters
        ----------
        hist : ROOT.TH2D or uproot.behaviors.TH2.TH2
        name : str
            name of the projected histogram
        axis : str
            `x` or `y`
        range : tuple
            range of the projection
        Returns
        -------
        ROOT.TH1D or uproot.behaviors.TH1.TH1
        Example
        -------
        # for ROOT.TH2D
        h2 = ROOT.TH2D('test', '', 100, -1., 1., 100, -1., 1.)
        for _ in range(100):
            h2.Fill(ROOT.gRandom.Gaus(0.), ROOT.gRandom.Gaus(0.5))
        h1 = hist_reader.projection(h2, name='h1', axis='x', range=(-0.5, 0.5))
        """
        axes = ['x', 'y']
        if axis.lower() not in axes:
            raise ValueError(f'axis must be either `x` or `y`, not {axis}')
        
        if name is None:
            name = hist.GetName() + f'_p{axis}'

        a = axes[1 - axes.index(axis.lower())]
        GetAxis = getattr(hist, f"Get{a.upper()}axis")
        xbin1 = GetAxis().FindBin(range[0])
        xbin2 = GetAxis().FindBin(range[1])

        return getattr(hist, f'Projection{axis.upper()}')(name, xbin1, xbin2)

    def df_to_hist1d(self, df, name, range, bins):
        """ convert a pd.DataFrame to ROOT.TH1D
        Parameters
        ----------
        df : pd.DataFrame
            must have columns `x`, `content`, `error`
        name : str
            name of the histogram
        range : tuple
            range of the histogram
        bins : int
            number of bins
        Returns
        ------- 
        ROOT.TH1D
        """
        h = ROOT.TH1D(name, '', bins, *range)
        for i in range(bins):
            h.SetBinContent(i+1, df['content'].iloc[i])
            h.SetBinError(i+1, df['error'].iloc[i])
        return h