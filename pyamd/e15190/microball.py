import pathlib
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from pyamd import DATABASE
from pyamd.utilities import ame
ame_table = ame.AME()

class microball:
    """ A class to handle the geometry and acceptance (in MeV for each ring) of the microball detector.
    """
    MAX_A = 100
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

        # a dictionary of dataframe to store the acceptance threshold for each ring
        self.threshold = dict()

    def configurate(self, path=None, reaction='Ca48Ni64E140'):
        """ Configuare the `ring`and `det` used in the experiment.
        Parameters
        ----------
        path : str
            Path to the configuration file. The file has the following format: `system`, `ring` and a list of `det` separated by space.
        reaction : str
            Reaction name, e.g. `Ca48Ni64E140`
        """
        if path is None:
            path = pathlib.Path(DATABASE, 'e15190/microball/acceptance/config.dat')

        self.configuration = dict()
        with open(str(path), 'r') as f:
            content = f.readlines()
            for line in content[1:]:
                if not reaction == line.split()[0]:
                    continue
                ring = line.split()[1]
                dets = list(map(int, line.split()[2:]))
                self.configuration[ring] = dets
        return self.configuration

    def _read_geometry(self, path, coordinate):
        """ Read the geometry file.
        Parameters
        ----------
        path : str
            Path to the geometry file, the file contains the following columns: `ring, det, theta_min, theta_max, phi_min and phi_max`
        coordinate : str
            'uball' : default
            'hira' : rotate :math: `\phi` angle by 90 degrees.  
        Returns
        -------
        df : pandas.DataFrame
            A pandas DataFrame contains the geometry information.
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

        for ring, dets in self.configuration.items():
            DET_MAX = 14
            for det in [det for det in range(1, DET_MAX + 1) if not det in dets]:
                if (ring, det) in df.index:
                    df.drop((ring, det))
        return df
    
    def get_theta_range(self, ring, det):
        return self.geometry.loc[ring, det]['theta_min'], self.geometry.loc[ring, det]['theta_max']

    def get_phi_range(self, ring, det):
        return self.geometry.loc[ring, det]['phi_min'], self.geometry.loc[ring, det]['phi_max']
    
    def get_ring(self, theta):
        """ Get the ring number from `self.geometry` given theta. """
        return self.geometry.query('theta_min <= @theta < theta_max').index.get_level_values('ring').values[0]
    def get_ring_and_detectorID(self, theta, phi):
        """ Get the ring number and detector ID from `self.geometry` given theta and phi. """
        return self.geometry.query('theta_min <= @theta < theta_max and phi_min <= @phi < phi_max').index.values[0]


    def set_threshold(self, ring=2, path=None, maxA=None):
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
        A = np.arange(1, maxA+1)
        # f = lambda x, a, b, c: a * x ** 2 + b * x + c
        f = lambda x, a, b: a * x + b
        
        threshold = dict()
        for Z, subdf in df.groupby('Z'):
            x = subdf['A'].to_numpy()
            y = subdf['kinergy_MeV'].to_numpy()

            popt, _ = curve_fit(f, x, y)
            missing_A = np.array([a for a in A if not a in x and ame_table.is_physical(a-Z, Z)])
            threshold[Z] = pd.concat([subdf, pd.DataFrame({
                'A': missing_A,
                'Z': Z,
                'kinergy_MeV': f(missing_A, *popt)
            })])
        self.threshold[ring] = pd.concat(threshold.values())

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
        return query['kinergy_MeV'].values[0]
    
    def is_inside(self, theta, phi):
        """ Check if the specified theta and phi are inside the microball. """
        return self.get_ring_and_detectorID(theta, phi) is not None

    def is_punchthrough(self, ring, A, Z, E):
        """ Check if the specified energy is detected by the specified ring. """
        return E >= self.get_threshold(ring, A, Z)

    def is_detected(self, ring, A, Z, E, theta, phi):
        """ Check if the particle is detected by microball. """
        return self.is_inside(theta, phi) and self.is_punchthrough(ring, A, Z, E)