import pyamd
from pyamd.analysis import centrality
from pyamd.utilities import root6, physics, style
from pyamd.e15190 import microball
import inspect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
style.set_matplotlib_style(mpl)

uball = microball.microball()


class MultiplicityAnalysis:

    def __init__(self, reaction='Ca48Ni64E140', skyrme='SkM', impact_parameter=(0., 10.), mode='all'):

        self.reaction = reaction
        self.skyrme = skyrme
        self.mode = mode
        self.impact_parameter = impact_parameter

        self.spectra = centrality.Multiplicity_ImpactParameter(
            reaction=reaction, skyrme=skyrme, impact_parameter=impact_parameter, mode=mode)

    def savefig(self, fig, name=''):
        name = inspect.cleandoc(f'''
            {name}_{self.reaction}_{self.skyrme}_{self.mode}_bmin{self.impact_parameter[0]}_bmax{self.impact_parameter[1]}.png
        ''')
        fig.savefig(name)

    def plot2d(self):
        fig, ax = plt.subplots()
        self.spectra.plot2d(ax, norm=LogNorm())
        ax.set(xlabel='Nc', ylabel='b [fm]')
        plt.tight_layout()
        self.savefig(fig, 'h2_multi_b')

    def plot_multiplity(self):
        fig, ax = plt.subplots()
        self.spectra.plot_multiplicity(ax, fmt='ko')
        ax.set(xlabel='Nc', yscale='log')
        plt.tight_layout()
        self.savefig(fig, 'h1_multi')

    def plot_impact_parameter(self, bins=20):
        fig, ax = plt.subplots()
        self.spectra.plot_impact_parameter(ax, bins=bins, fmt='ko')
        ax.set(xlabel='b [fm]', ylabel='normalized dM/db [fm-1]', yscale='log')
        plt.tight_layout()
        self.savefig(fig, 'h1_impact_parameter')

    def plot_dsigma_db(self):
        fig, ax = plt.subplots()
        self.spectra.plot_dsigma_db(ax, unit='mb', fmt='ko')
        self.spectra.plot_fitted_dsigma_db(ax, unit='mb', fmt='r--')
        ax.set(xlabel='b [fm]', ylabel='dsigma/db [mb fm-1]')
        plt.tight_layout()
        self.savefig(fig, 'h1_dsigma_db')

    def plot_multiplicity_slices(self):
        colors = plt.rcParams["axes.prop_cycle"]()
        ranges = np.array([0., 1.5, 2.5, 3., 4., 5., 6., 7., 8., 9.])
        ranges = list(zip(ranges[:-1], ranges[1:]-0.01))
        fig, axes = plt.subplots(3, 3, figsize=(
            8, 6), sharex=True, sharey=True)
        for range, ax in zip(ranges, axes.flatten()):

            c = next(colors)["color"]
            self.spectra.plot_multiplicity(
                ax, cut=range, label=f'b=[{range[0]:.1f}-{range[1]:.1f})', color=c)
            self.spectra.plot_fitted_multiplicity(
                ax, cut=range, label='gaus. fit', color=c, fmt='--')
            ax.set(xlim=(-0.5, 25.5), ylim=(1e-6, 5e-1), yscale='log')
            ax.legend()
        fig.supxlabel('Nc')
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.tight_layout()
        self.savefig(fig, 'h1_multiplicity_different_bimp')

    def plot_impact_parameter_slices(self):
        colors = plt.rcParams["axes.prop_cycle"]()
        ranges = np.array([1, 2, 3, 4, 6, 8, 10, 12, 14, 25])
        ranges = list(zip(ranges[:-1], ranges[1:]-0.01))
        fig, axes = plt.subplots(3, 3, figsize=(
            8, 6), sharex=True, sharey=True)
        for range, ax in zip(ranges, axes.flatten()):
            c = next(colors)["color"]
            self.spectra.plot_impact_parameter(
                ax, cut=range, label=f'Nc=[{range[0]:.0f}-{range[1]:.0f})', color=c, fmt='.')
            self.spectra.plot_fitted_impact_parameter(
                ax, cut=range, label='gaus. fit', color=c, fmt='--')
            ax.set(xlim=(0., 10.), ylim=(1e-6, 5e-1), yscale='log')
            ax.legend()
        fig.supxlabel('b [fm]')
        plt.tight_layout()
        self.savefig(fig, 'h1_impact_parameter_different_multiplicity')

    def plot_mapping(self):
        fig, ax = plt.subplots()
        self.spectra.plot_mapping(ax, fmt='r--', label='AMD')
        uball.plot_impact_parameter_mapping(
            ax, path=f'{str(pyamd.PROJECT_DIR)}/database/microball/bimp_mapping/{self.reaction}.dat', fmt='k.', label="e15190")
        ax.set(xlabel='Nc [fm]', ylabel='b [fm]')
        ax.legend()
        plt.tight_layout()
        self.savefig(fig, 'h1_mapping')

    def run(self):

        self.plot2d()
        self.plot_multiplity()
        self.plot_impact_parameter()
        self.plot_dsigma_db()
        self.plot_multiplicity_slices()
        self.plot_impact_parameter_slices()
        self.plot_mapping()


if __name__ == '__main__':
    anal = MultiplicityAnalysis(
        reaction='Ca40Ni58E140', skyrme='SkM', impact_parameter=(0., 10.), mode='all')
    anal.run()
