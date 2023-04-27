from pyamd import PROJECT_DIR

import pandas as pd
import pathlib
import requests
import re
import fortranformat as ff
from astropy import units, constants

path_download = pathlib.Path(PROJECT_DIR, 'database/ame')


class AME:
    def __init__(self, url=None, dl_path=None):
        if url is None:
            url = 'https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt'
        r = requests.get(url)

        if dl_path is None:
            path_download.mkdir(exist_ok=True, parents=True)
            dl_path = pathlib.Path(path_download, 'mass20.txt')
        self.dl_path = dl_path

        with open(str(dl_path), 'w') as file:
            file.write(r.text)

        self._read()


    def _read(self):
        with open(str(self.dl_path), 'r') as file:
            lines = file.readlines()

        for line in lines:
            if 'format' in line:
                format = line.split(':')[1].strip()
                break
        for i, line in enumerate(lines):
            if 'keV' in line:
                break
        lines = lines[i+1:]
        format = '(' + format + ')'

        freader = ff.FortranRecordReader(format)

        columns_names = [
            '1',
            'N-Z',
            'N',
            'Z',
            'A',
            'EL',
            'O',
            'mass_excess',
            'mass_excess_err',
            'binding_energy_per_nucleon',
            'binding_energy_per_nucleon_err',
            'B-',
            'beta_decay_energy',
            'beta_decay_energy_err',
            'au',
            'atomic_mass_au',
            'atomic_mass_err',
        ]

        # '#' means estimated value (non-experimental)
        # '*' means unable to calculate

        # fortranformat does not handle '*'
        lines = [re.sub('[*]', ' ', line) for line in lines]  

        # avoid unit change in fortranformat
        lines = [re.sub('[#]', '.', line) for line in lines]  

        lines = [freader.read(line) for line in lines]
        
        df = pd.DataFrame(lines, columns=columns_names)

        df['symbol'] = list(map(lambda s: s.lower().strip(), df['EL']))
        df['symbol'] = df['symbol'] + df['A'].astype(str)

        df['N'] = df['N'].astype(int)
        df['Z'] = df['Z'].astype(int)
        df['A'] = df['A'].astype(int)

        # change of units

        df['mass_excess'] = df['mass_excess'] * units.keV.to('MeV', 1.)
        df['mass_excess_err'] = df['mass_excess_err'] * units.keV.to('MeV', 1.)

        df['au'] = df['A'] + df['atomic_mass_au'] * units.micron.to('m', 1.)
        df['au_err'] = df['atomic_mass_err'] * units.micron.to('m', 1.)

        df['mass'] = df['au'] * (units.u * constants.c ** 2).to('MeV')
        df['mass_err'] = df['au_err'] * (units.u * constants.c ** 2).to('MeV')

        df['binding_energy_per_nucleon'] = df['binding_energy_per_nucleon'] * units.keV.to('MeV', 1.)
        df['binding_energy_per_nucleon_err'] = df['binding_energy_per_nucleon_err'] * units.keV.to('MeV', 1.)

        df.drop(columns=['1', 'N-Z', 'O', 'B-', 'EL', 'atomic_mass_au', 'atomic_mass_err'], inplace=True)

        self.df = df

    def _get_row(self, symbol=None, N=None, Z=None):
        if N is not None and Z is not None:
            return self._get_row_from_NZ(N, Z)
        elif symbol is not None:
            return self._get_row_from_symbol(symbol)
        else:
            raise ValueError('Either symbol or N and Z must be specified.')

    def _get_row_from_symbol(self, symbol):
        return self.df.loc[self.df['symbol'] == symbol.lower()]
    
    def _get_row_from_NZ(self, N, Z):
        return self.df.loc[((self.df['N'] == N) & (self.df['Z'] == Z))]

    def get_symbol(self, N, Z):
        return self._get_row(N=N, Z=Z)['symbol'].values[0]
    
    def get_NZ(self, symbol):
        row = self._get_row(symbol)
        return (row['N'].values[0], row['Z'].values[0])

    def get_mass(self, symbol='', N=None, Z=None, unit='MeV'):
        row = self._get_row(symbol, N, Z)
        return row['mass'].values[0] * units.MeV.to(unit)
    
    def get_binding_energy(self, symbol=None, N=None, Z=None, unit='MeV'):
        return self._get_row(symbol, N, Z)['binding_energy_per_nucleon'].values[0] * units.MeV.to(unit)
    
    def get_mass_data(self, fname='ame_mass.txt', usecols=['symbol', 'Z', 'A', 'mass', 'binding_energy_per_nucleon']):
        """ Output the mass table to a format convenient for reading in C++.
        Parameters
        ----------
        fname : str
            Output file name.
        usecols : list
            Columns to be output.
        """
        if not set(usecols).issubset(self.df.columns):
            raise ValueError('usecols must be a subset of the columns of the AME table.')
        
        df = self.df[usecols].copy()
        df.to_csv(fname, sep=' ', columns=usecols, index=False)

        return