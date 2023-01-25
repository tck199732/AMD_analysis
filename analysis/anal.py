from pyamd import PROJECT_DIR, DATABASE

import subprocess
import pickle
import pathlib
import time

"""
This script works for data obtained in Dec2021. The names of required primary files are saved in `amdfiles`, the pickled dictionary.
"""

reaction = 'Ca40Ni58E56'
skyrme = 'SkM'
Nc = 10

amdfiles = pathlib.Path(DATABASE, 'data/amdfiles_dec2021.pkl')
input_dir = pathlib.Path(DATABASE, 'data/amd/dec2021/b3fm/filtered')
out_dir = pathlib.Path(PROJECT_DIR, f'result/spectra/sigma0.85_b3fm')

def main():

    out_dir.mkdir(exist_ok=True, parents=True)
    
    exe = pathlib.Path(PROJECT_DIR, 'analysis/anal.exe')
    if not exe.exists():
        subprocess.run('make', shell=True)
    
    path_filtered_table3 = pathlib.Path(
        input_dir, f'{reaction}_{skyrme}_table3.root')

    with open(str(amdfiles), 'rb') as f:
        primary_files = pickle.load(f)
    
    path_raw_table21 = ' '.join([str(pth.with_suffix('.root')) for pth in primary_files[reaction][skyrme]['primary']])

    path_out = pathlib.Path(
        out_dir, f'{reaction}_{skyrme}_Nc{Nc}.root')

    args = ' '.join([
        '-r' , reaction,
        '-o' , str(path_out), 
        '-f' , str(path_filtered_table3),
        '-p' , f'\"{path_raw_table21}\"',
        '-c' , str(Nc),
    ])
    
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')


if __name__ == '__main__':
    main()
