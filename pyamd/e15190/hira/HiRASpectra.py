

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyamd.analysis import spectra
from pyamd.utilities import helper, minuit
from pyamd.e15190 import e15190
from pyamd import PROJECT_DIR
import pathlib
import functools
import pandas as pd
import numpy as np
import random
from datetime import datetime
random.seed(datetime.now().timestamp())


class HiRASpectra:
    DIR = pathlib.Path(PROJECT_DIR, 'database/hira/spectra')

    def __init__(self, reaction, particle, path=None):
        self.reaction = reaction
        reaction = e15190.reaction(reaction)
        self.betacms = reaction.get_betacms()
        self.beam_rapidity = reaction.get_rapidity_beam()

        self.particle = particle
        if path is None:
            path = self.DIR / self.reaction / f'h2_pta_rapidity_{particle}.dat'

        df = pd.read_csv(str(path), delim_whitespace=True)
        self.PtRapidityLAB = spectra.PtRapidityLAB(
            particle=self.particle, df=df, reaction=self.reaction)

    def PtSpectrum(self, xrange=(0.4, 0.6), yrange=(0, 600), bins=30, correct_coverage=False, correct_range=(0, 600)):
        return self.PtRapidityLAB.PtSpectrum(xrange, yrange, bins, correct_coverage, correct_range)

    def plotPtSpectrum(self, ax=None, xrange=(0.4, 0.6), yrange=(0, 600), bins=30, correct_coverage=False, correct_range=(0, 600), **kwargs):
        return self.PtRapidityLAB.plotPtSpectrum(ax=ax, xrange=xrange, yrange=yrange, bins=bins, correct_coverage=correct_coverage, correct_range=correct_range, **kwargs)

    def plotRapidityCMS(self, ax=None, range=(-0.5, 0.5), bins=100, **kwargs):
        return self.PtRapidityLAB.plotRapidityCMS(ax=ax, range=range, bins=bins, **kwargs)

    def plotPtRapidity(self, ax=None, range=[[0., 1.], [0., 600.]], bins=[100, 600], **kwargs):
        return self.PtRapidityLAB.plotPtRapidity(ax=ax, range=range, bins=bins, **kwargs)

    # @functools.cached_property
    def SetEkinThetaCMS(self, xrange=(0.4, 0.6), yrange=(0., 600.), correct_coverage=False):
        df = self.PtRapidityLAB.EkinThetaCMS(
            xrange=xrange, yrange=yrange, correct_coverage=correct_coverage)
        self.EkinThetaCMS = spectra.EkinThetaCMS(
            particle=self.particle, df=df)

    def EkinCMS(self, xrange=(0., 200.), yrange=(0., 180.), bins=50):
        return self.EkinThetaCMS.EkinCMS(xrange, yrange, bins)

    def plotEkinCMS(self, ax=None, range=(0, 200), bins=50, **kwargs):
        return self.EkinThetaCMS.plotEkinCMS(ax=ax, range=range, bins=bins, **kwargs)

    def plotThetaCMS(self, ax=None, range=(40, 140), bins=100, **kwargs):
        return self.EkinThetaCMS.plotThetaCMS(ax=ax, range=range, bins=bins, **kwargs)

    def plotEkinThetaCMS(self, ax=None, range=[[0, 200], [40, 140]], bins=[50, 100], **kwargs):
        return self.EkinThetaCMS.plotEkinThetaCMS(ax=ax, range=range, bins=bins, **kwargs)


if __name__ == '__main__':

    hira = HiRASpectra('Ca48Ni64E140', 'p')

    fig, ax = plt.subplots()
    hira.plotPtRapidity(ax, norm=LogNorm())
    fig.savefig('h2_pta_rapidity_p.png')

    fig, ax = plt.subplots()
    hira.plotPtSpectrum(ax)
    fig.savefig('h1_pta_p.png')

    fig, ax = plt.subplots()
    hira.plotRapidityCMS(ax)
    fig.savefig('h1_rap_p.png')
