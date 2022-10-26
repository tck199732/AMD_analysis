
import pathlib
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

try:
    import ROOT
except:
    path_root = subprocess.check_output('which root', shell=True)
    path_root = str(path_root, encoding='utf-8').strip()
    path_root = pathlib.Path(str(path_root)).parent.parent
    path_root = pathlib.Path(path_root, 'lib')
    sys.path.append(str(path_root))
    import ROOT


class TFile:
    def __init__(self, path, mode='READ'):
        self.path = pathlib.Path(path)
        self.mode = mode.upper()

    def __enter__(self):
        self.file = ROOT.TFile(str(self.path), self.mode)
        return self.file

    def __exit__(self, *args):
        self.file.Close()


class HistogramReader:
    def __init__(self):
        pass

    def hist1d_to_df(self, hist, xname='x', yname='y', keep_zeros=True):
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

        return df if keep_zeros else df.query(f'{zname} != 0.0').reset_index(drop=True)

    def hist2d_to_df(self, histo, xname='x', yname='y', zname='z', keep_zeros=True):
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

    def graph_to_df(self, grer, xname='x', yname='y', keep_zeros=True):
        df = dict()
        columns = [xname, yname, f'{yname}_err']
        getters = [grer.GetX, grer.GetY, grer.GetEY]
        for col, get in zip(columns, getters):
            df[col] = [get()[b] for b in range(grer.GetN())]

        df = pd.DataFrame(df)
        condition = (df[yname] != 0.0)
        df[f'{yname}_ferr'] = np.where(
            condition, np.abs(df[f'{yname}_err']/df[yname]), 0.0)

        return df if keep_zeros else df.query(f'{zname} != 0.0').reset_index(drop=True)
