import pathlib
import numpy as np
import collections
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from typing import Literal

from pyamd import DATABASE
from pyamd.utilities import style
style.set_matplotlib_style(mpl)


class microball:

    def __init__(self, path=None, coordinate='uball'):
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/geometry.dat')
        
        self.coordinate = coordinate
        self.geometry = self._read_geometry(path, coordinate)
        self.configuration = dict()
        self.threshold = collections.defaultdict(dict)

        self.configurated = False
    
    def _read_geometry(self, path, coordinate):
        """ Base geometry of microball.  
        Parameter
        ---------
        coordinate : str
            'uball' : default
            'hira' : rotate :math: `\phi` angle by 90 degrees.
        """
        df = pd.read_csv(str(path), delim_whitespace=True).set_index(['ring', 'det'])
        if coordinate.lower() == 'hira':
            df['phi_min'] += 90.
            df['phi_max'] += 90.
            
        for i, row in df.iterrows():
            if row['phi_min'] >= 360:
                df.at[i, 'phi_min'] -= 360

            if row['phi_max'] > 360:
                df.at[i, 'phi_max'] -= 360

        return df
    
    def get_theta_range(self, ring, det):
        return self.geometry.loc[ring, det]['theta_min'], self.geometry.loc[ring, det]['theta_max']

    def get_phi_range(self, ring, det):
        return self.geometry.loc[ring, det]['phi_min'], self.geometry.loc[ring, det]['phi_max']

    
    def configurate(self, path=None, reaction='Ca48Ni64E140'):
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/config.dat')
        with open(str(path), 'r') as f:
            content = f.readlines()
            for line in content[1:]:
                if not reaction == line.split()[0]:
                    continue

                ring = line.split()[1]
                dets = list(map(int, line.split()[2:]))
                self.configuration[ring] = dets
                DET_MAX = 14
                for det in [det for det in range(1, DET_MAX + 1) if not det in dets]:
                    if (ring, det) in self.geometry.index:
                        self.geometry.drop((ring, det))
        
        self.configurated = True
        return self.geometry

                
    def plot_coverage(self, cmap='viridis', **kwargs):
        
        figdict = dict(
            figsize = (6,4),
        )
        figdict.update(kwargs)
        
        fig, ax = plt.subplots(**figdict)
        cm = plt.get_cmap(cmap)(np.linspace(0.3, 1, len(self.geometry)))

        for ir, (index, row) in enumerate(self.geometry.iterrows()):
            ring, det = index
            x = row['phi_min']
            y = row['theta_min']
            dx = row['phi_max'] - x 
            dy = row['theta_max'] - y
            patch = plt.Rectangle(
                (x, y), dx, dy, 
                fill=True, alpha=0.6, edgecolor='k',
                facecolor=cm[ir], linewidth=1
            )
            ax.add_patch(patch)
            ax.text(
                x + 0.5 * dx, 
                y + 0.5 * dy, f'{ring}:{det}',
                ha='center', va='center', 
                fontdict={
                    'size': 8, 
                    'weight': 'roman', 
                    'family' : 'serif', 
                    'style' : 'normal',
                } 
            )

        ax.set(title=r'E15190 $\mu$-ball coverage')
        ax.set(xlabel=r'$\phi$ (deg.)', ylabel=r'$\theta$ (deg.)')
        if self.coordinate == 'uball':
            ax.set(xlim=(-30,360), ylim=(0,160))
        elif self.coordinate == 'hira':
            ax.set(xlim=(0,360), ylim=(0,160))

        plt.tight_layout()

        return fig, ax


    def set_threshold(self, ring=2, path=None):
        """
        set_threshold(2, threshold_Sn_65mgcm2.dat)
        set_threshold(3, threshold_Sn_58mgcm2.dat)
        set_threshold(4, threshold_Sn_50mgcm2.dat)
        set_threshold(5, threshold_Sn_43mgcm2.dat)
        set_threshold(7, threshold_Sn_30mgcm2.dat)
        set_threshold(8, threshold_Sn_23mgcm2.dat)
        """
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/threshold_Sn_65mgcm2.dat')
        
        df = pd.read_csv(str(path), delim_whitespace=True)
        df.columns = ['A', 'Z', 'threshold']
        df.set_index(['A', 'Z'])
        
        for AZ, row in df.iterrows():
            self.threshold[ring] = float(row[AZ]['threshold'])
        
        # here, do a fit to extend the data, to be done later.

        return

    def get_threshold(self, ring, A, Z):
        
        if not ring in self.threshold:
            raise ValueError('Not yet set up threshold data.')
        
        return self.threshold[ring][(A,Z)]


class MultiplicityMapping:
    def __init__(self, path=None, reaction='Ca48Ni64E140'):
        if path is None:
            path = pathlib.Path(DATABASE, f'e15190/microball/bimp_mapping/{reaction}.dat')
        
        self.df = pd.read_csv(str(path), delim_whitespace=True).reset_index(drop=True)
    
    def MultiplicityToImpactParameter(self, m=5, y:Literal['b', 'bhat']='b'):
        df = self.df.loc[self.df['multiplicity']==m]
        return (df[y].values[0], df[f'{y}_err'].values[0])

