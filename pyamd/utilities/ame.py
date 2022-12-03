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

        self.read()

    def read(self):
        with open(str(self.dl_path), 'r') as file:
            lines = file.readlines()

        capture_format = False
        counter_0 = 0
        for i, line in enumerate(lines):
            if 'format' in line:
                if capture_format is False:
                    format = line.split(':')[1].strip()
                    capture_format = True
            if line.startswith('0'):
                counter_0 += 1
            if counter_0 == 3:
                break

        lines = lines[i:]

        self.format = f'({format})'
        freader = ff.FortranRecordReader(self.format)

        columns_names = {
            2: 'N',
            3: 'Z',
            4: 'A',
            5: 'symbol',
            7: 'mass_excess',
            8: 'mass_excess_err',
            15: 'atomic_mass',
            16: 'atomic_mass_err',
        }

        # '#' means estimated value (non-experimental)
        # '*' means unable to calculate
        lines = [re.sub('[*]', ' ', line)
                 for line in lines]  # fortranformat does not handle '*'
        lines = [re.sub('[#]', '.', line)
                 for line in lines]  # avoid unit change in fortranformat
        lines = [freader.read(line) for line in lines]
        for i, line in enumerate(lines):
            lines[i] = [line[id] for id in columns_names]
            # print(lines[i])

        df = pd.DataFrame(lines, columns=columns_names.values())
        df['N'] = df['N'].astype(int)
        df['Z'] = df['Z'].astype(int)
        df['A'] = df['A'].astype(int)

        df['symbol'] = list(map(lambda s: s.lower().strip(), df['symbol']))
        df['symbol'] = df['symbol'] + df['A'].astype(str)
        # df['symbol'] = list(map(lambda s: s.lower*(), df['symbol']))

        df['atomic_mass'] = df['A'] + \
            df['atomic_mass'] * units.micron.to('m', 1.)
        df['atomic_mass_err'] = df['atomic_mass_err'] * \
            units.micron.to('m', 1.)

        # df['mass'] = df['A'] + df['atomic_mass'] * units.micron.to(units.m, 1.)
        # df['mass'] = df['mass'] * constants.u
        # df['mass'] = list(
        # map(lambda m: (m * units.kg * constants.c**2).to('MeV').value, df['mass']))
        # print(df)
        self.df = df

    def get_NZ(self, symbol):
        row = self.df.loc[self.df['symbol'] == symbol.lower()]
        if row.empty:
            raise ValueError(f'symbol not found : {symbol}')
        return (row['N'].values[0], row['Z'].values[0])

    def get_mass(self, symbol='', N=None, Z=None):
        if N is None or Z is None:
            N, Z = self.get_NZ(symbol)
        row = self.df.loc[((self.df['N'] == N) & (self.df['Z'] == Z))]
        if row.empty:
            if symbol is None:
                raise ValueError(f'symbol not found : {symbol}')
            elif N is None or Z is None:
                raise ValueError('please provide N and Z.')

        # atomic_mass = row['atomic_mass'].values[0]
        # return (atomic_mass * constants.u * constants.c**2).to('MeV')
        mass_excess = row['mass_excess'].values[0]
        A = row['A'].values[0]
        return (A * constants.u * constants.c**2 + mass_excess * units.keV).to('MeV')

    def get_symbol(self, N, Z):
        row = self.df.loc[((self.df['N'] == N) & (self.df['Z'] == Z))]
        return row['symbol'].values[0]


if __name__ == '__main__':
    ame = AME()
    print(ame.get_NZ('ca40'))
    print(ame.get_mass('h1'))
