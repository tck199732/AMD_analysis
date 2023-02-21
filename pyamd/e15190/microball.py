import pathlib
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from pyamd import DATABASE
from pyamd.utilities import style
style.set_matplotlib_style(mpl)

class microball:
    """ A class to handle the geometry and acceptance (in MeV for each ring) of the microball detector.
    """
    MAX_A = 100
    MAX_Z = 100
    def __init__(self, path=None, coordinate='uball'):
        """
        Parameters
        ----------
        path : str
            Path to the geometry file, the file contains the following columns: ring, det, theta_min, theta_max, phi_min and phi_max
        coordinate : str
            'uball' : default
            'hira' : rotate :math: `\phi` angle by 90 degrees.
        """
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/geometry.dat')
        
        self.coordinate = coordinate
        self.geometry = self._read_geometry(path, coordinate)
        self.configuration = dict()
        self.threshold = dict()

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
        """ configurate according to the `ring` and `det` used in the experiment.
        Parameters
        ----------
        path : str
            Path to the configuration file. The file has the following format: `system`, `ring` and a list of `det` separated by space.
        reaction : str
            Reaction name, e.g. `Ca48Ni64E140`
        """
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

    def set_threshold(self, ring=2, path=None, maxA=None, maxZ=None):
        """ Set the threshold energy (MeV) for the specified ring according to the threshold file. In the file, the first column is the `A` of the target, the second column is the `Z` of the target, and the third column is the threshold energy (MeV). If the `A` and `Z` of the target are not found in the file, the threshold energy will be interpolated from the nearest `A` and `Z` values. 
        Parameters
        ----------
        ring : int
            Ring number.
        path : str
            Path to the threshold file. The file has the following format: `A`, `Z` and `threshold` separated by space.
        maxA : int, optional
            Maximum `A` of the target. 
        maxZ : int, optional
            Maximum `Z` of the target.
        Examples
        --------
        set_threshold(2, threshold_Sn_65mgcm2.dat)
        set_threshold(3, threshold_Sn_58mgcm2.dat)
        set_threshold(4, threshold_Sn_50mgcm2.dat)
        set_threshold(5, threshold_Sn_43mgcm2.dat)
        set_threshold(7, threshold_Sn_30mgcm2.dat)
        set_threshold(8, threshold_Sn_23mgcm2.dat)
        """
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/threshold_Sn_65mgcm2.dat')
        
        df = pd.read_csv(str(path), delim_whitespace=True, usecols=[0,1,2])
        df.columns = ['A', 'Z', 'kinergy_MeV']

        maxA = self.MAX_A if maxA is None else maxA
        maxZ = self.MAX_Z if maxZ is None else maxZ

        for Z, subdf in df.groupby('Z'):
            a = subdf['A'].to_numpy()
            e = subdf['kinergy_MeV'].to_numpy()
            # interpolate and extrapolate the threshold energy
        

        # setting self.thereshld[ring] = df

    def get_threshold(self, ring, A, Z):
        if not ring in self.threshold:
            raise ValueError('Not yet set up threshold data.')
        if not A >= Z:
            raise ValueError('A must be greater than or equal to Z.')
        if A < 0 or Z < 0:
            raise ValueError('A and Z must be greater than or equal to 0.')
        
        query = self.threshold[ring].query('A == @A and Z == @Z')
        if len(query) == 0:
            raise ValueError('Not found the specified target.')
        if len(query) > 1:
            warnings.warn('Found more than one target. Use the first one.')
        return query['threshold'].values[0]
    
if __name__ == '__main__':
    mb = microball()
    mb.configurate()
    mb.set_threshold(2, path=pathlib.Path(DATABASE, 'e15190/microball/acceptance/threshold_Sn_65mgcm2.dat'))