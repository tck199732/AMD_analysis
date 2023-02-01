import pickle
import glob
import pathlib
import itertools
import re
from pyamd import PROJECT_DIR, DATABASE
from pyamd.utilities import misc
"""
This script collects AMD files as a nested dictionary. The naming convention of AMD files is `{beam}{target}_En{energy}MeV_{skyrme}_{statistics}k_{table}x.dat`, e.g. `Ca48Ni64_En140MeV_SkM_110k_table3.dat`. User provides the input directory containing these files, this program will collect the file names sorted by the statistics. For instance, ['Ca48Ni64E140']['SkM']['secondary'] = ['/data/Ca48Ni64_En140MeV_SkM_110k_table3.dat', '/data/Ca48Ni64_En140MeV_SkM_150k_table3.dat','/data/Ca48Ni64_En140MeV_SkM_300k_table3.dat']
"""

output_path = pathlib.Path(DATABASE, 'data/amdfiles_dec2021.pkl')
dat_dir = pathlib.Path(PROJECT_DIR, 'database/data/amd/dec2021/b3fm')

##########################################################################
beam = ['Ca40', 'Ca48']
target = ['Ni58', 'Ni64', 'Sn112', 'Sn124']
energy = [56, 140]
skyrme = ['SkM', 'SLy4', 'SLy4_L108']
##########################################################################


def main():

    files_list = misc.InfiniteDict()
    for b, t, e, sky in itertools.product(beam, target, energy, skyrme):
        if (beam.index(b) + target.index(t)) % 2 != 0:
            continue

        for operation_mode in ['table21', 'table3']:
            fnames = [pathlib.Path(f).stem for f in glob.iglob(
                f'{str(dat_dir)}/{b}{t}_En{e}MeV_{sky}_*_{operation_mode}.dat')]

            # remove cases like SLy4_L108 from SLy4
            fnames = [f'{f}.root' for f in fnames if re.match(
                '[0-9]+k', f.removeprefix(f'{b}{t}_En{e}MeV_{sky}_').split('_')[0])]

            # sort according to statistics
            stats = set([re.findall('[0-9]+k', n)[0] for n in fnames])
            stats = sorted([int(s[:-1]) for s in stats])

            fnames = [pathlib.Path(dat_dir, f'{b}{t}_En{e}MeV_{sky}_{s}k_{operation_mode}.dat') for s in stats]

            # assign key name to table21 and table3
            if operation_mode == 'table21':
                key = 'primary'
            elif operation_mode == 'table3':
                key = 'secondary'

            files_list[f'{b}{t}E{e}'][sky][key] = fnames

    with open(str(output_path), 'wb') as f:
        pickle.dump(files_list, f)
    

if __name__ == '__main__':
    main()

    # testing
    with open(str(output_path), 'rb') as f:
        data = pickle.load(f)
        print(data['Ca48Ni64E140']['SLy4']['primary'])