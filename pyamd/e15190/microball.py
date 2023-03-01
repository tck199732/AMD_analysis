import pathlib
import warnings
import inspect
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

        for ring, dets in self.configuration.items():
            DET_MAX = 14
            for det in [det for det in range(1, DET_MAX + 1) if not det in dets]:
                if (ring, det) in self.geometry.index:
                    self.geometry.drop((ring, det))

        return
    
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

    def is_detected(self, A, Z, E, theta, phi):
        """ Check if the particle is detected by microball. """
        ring = self.get_ring(theta)
        return self.is_inside(theta, phi) and self.is_punchthrough(ring, A, Z, E)
    
    def is_inside_cpp(self, theta_br='theta_deg', phi_br='phi_deg', fcn_name='is_inside'):
        """ Construct a cpp function code snipet to be fed to `ROOT.RDataFrame.Define`
        """
        conditions = ' ||\n'.join([inspect.cleandoc(f'''
        (theta >= {row['theta_min']} && theta < {row['theta_max']} && phi >= {row['phi_min']} && phi < {row['phi_max']})''') for _, row in self.geometry.iterrows()])

        script = inspect.cleandoc(f''' 
            ROOT::RVec<bool> {fcn_name}(const ROOT::RVec<double>& {theta_br}, const ROOT::RVec<double>& {phi_br}) {{
                ROOT::RVec<bool> result;
                int n = {theta_br}.size();

                for (int i = 0; i < n; i++) {{
                    double theta = {theta_br}[i];
                    double phi = {phi_br}[i];
                    result.push_back({conditions});
                }}
                return result;
            }}
        ''')
        return script
    
    def get_ring_cpp(self, theta_br='theta_deg', fcn_name='get_ring'):

        df = self.geometry.reset_index()[['ring', 'theta_min', 'theta_max']].drop_duplicates().reset_index(drop=True)

        mapping = {
            row['ring'] : (row['theta_min'], row['theta_max']) for _, row in df.iterrows()
        }
        
        if_else_statement = ''
        for i, (key, value) in enumerate(mapping.items()):
            if i == 0:
                if_else_statement = f'''
                    if (theta >= {value[0]} && theta < {value[1]}) {{
                        ring = {key};
                    }}
                '''
            else:
                if_else_statement += f'''
                    else if (theta >= {value[0]} && theta < {value[1]}) {{
                        ring = {key};
                    }}
                '''
        script = f'''
            ROOT::RVec<int> {fcn_name}(const ROOT::RVec<double>& {theta_br}) {{
                ROOT::RVec<int> result;
                int n = {theta_br}.size();
                for (int i = 0; i < n; i++) {{
                    int ring = -1;
                    double theta = {theta_br}[i];
                    {if_else_statement}
                    result.push_back(ring);
                }}
                return result;
            }}
        '''
        return script
    

    def get_CsI_cpp(self, theta_br='thete_deg', phi_br='phi_deg', fcn_name='get_CsI'):
        df = self.geometry.reset_index()
        detector_id = [f"""{{ 
           {{ {row['theta_min']}, {row['theta_max']}, {row['phi_min']}, {row['phi_max']} }}, {{ {row['ring']:.0f}, {row['det']:.0f} }}
        }}""" for _, row in df.iterrows()]
        
        detector_id = ',\n'.join(list(map(inspect.cleandoc, detector_id)))
        

        script = f'''
            std::map<std::tuple<double, double, double, double>, std::tuple<int, int>> detector_id = {{
                {detector_id}
            }};

            ROOT::RVec<int> {fcn_name}(const ROOT::RVec<double>& {theta_br}, const ROOT::RVec<double>& {phi_br}) {{
                ROOT::RVec<int> result;
                int n = {theta_br}.size();
                for (int i = 0; i < n; i++) {{
                    int ring = -1;
                    int det = -1;
                    double theta = {theta_br}[i];
                    double phi = {phi_br}[i];
                    for (auto& [key, value] : detector_id) {{
                        auto [theta_min, theta_max, phi_min, phi_max] = key;
                        if (theta >= theta_min && theta < theta_max && phi >= phi_min && phi < phi_max) {{
                            std::tie(ring, det) = value;
                            break;
                        }}
                    }}
                    result.push_back(det);
                }}  
                return result;
            }}
        '''
        return script


    def get_threshold_data_cpp(self):

        df = self.threshold
        for r, subdf in df.items():
            subdf['ring'] = r
        df = pd.concat([subdf for subdf in self.threshold.values()])

        # construct a cpp map for the threshold data
        return ',\n'.join([inspect.cleandoc(f"""
            {{ {{ {row['ring']:.0f}, {row['A']:.0f}, {row['Z']:.0f}  }}, {row['kinergy_MeV']} }}
        """) for _, row in df.iterrows()])

    def generate_threshold_data(self, fname='threshold_data.dat'):
        df = self.threshold
        for r, subdf in df.items():
            subdf['ring'] = r
        df = pd.concat([subdf for subdf in self.threshold.values()])
        
        df.to_csv(fname, columns=['ring', 'A', 'Z', 'kinergy_MeV'], index=False, sep=' ')
        return

    def is_punchthrough_cpp(self, ring_br='ring', A_br='A', Z_br='Z', E_br='kinergy', fcn_name='is_punchthrough'):
        threshold_data = self.get_threshold_data_cpp()
        script = f'''
            std::map<std::tuple<int, int, int>, double> threshold_map = {{
                {threshold_data}
            }};
            ROOT::RVec<bool> {fcn_name}(const ROOT::RVec<int>& {ring_br}, const ROOT::RVec<int>& {A_br}, const ROOT::RVec<int>& {Z_br}, const ROOT::RVec<double>& {E_br}) {{
                ROOT::RVec<bool> result;
                int n = {ring_br}.size();
                for (int i = 0; i < n; i++) {{
                    int r = {ring_br}[i];
                    int a = {A_br}[i];
                    int z = {Z_br}[i];
                    double e = {E_br}[i];
                    
                    bool is_punchthr = threshold_map.count({{r, a, z}}) > 0  && e >= threshold_map[{{r, a, z}}];
                    result.push_back(is_punchthr);
                }}
                return result;
            }}
        '''
        return script

if __name__ == '__main__':
    mb = microball()
    mb.configurate(reaction='Ca48Ni64E140')
    # setting up the threshold data for microball
    threshold_map = {
        2: 'threshold_Sn_65mgcm2.dat',
        3: 'threshold_Sn_58mgcm2.dat',
        4: 'threshold_Sn_50mgcm2.dat',
        5: 'threshold_Sn_43mgcm2.dat',
        7: 'threshold_Sn_30mgcm2.dat',
        8: 'threshold_Sn_23mgcm2.dat',
    }
    for ring, filename in threshold_map.items():
        mb.set_threshold(ring, pathlib.Path(DATABASE, 'e15190/microball/acceptance') / filename)

    mb.generate_threshold_data()