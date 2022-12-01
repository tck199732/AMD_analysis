from pyamd import PROJECT_DIR, DATABASE, EXE_DIR, SRC_DIR
import subprocess
import pathlib
import os
import time

path_main = pathlib.Path(EXE_DIR, 'spectra.cpp')

# input_dir = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
# out_dir = pathlib.Path(PROJECT_DIR, 'result/spectra/sigma0.85_b3fm')

input_dir = pathlib.Path('/data/amd/nov2022/sigma0.85/filtered')
out_dir = pathlib.Path(PROJECT_DIR, 'result/spectra/sigma0.85_b3fm')
out_dir.mkdir(exist_ok=True, parents=True)


# reaction = 'Ca40Ni58'
reaction = 'Ca48Sn124'
energy = 140
skyrme = 'SkM'
Nc = (1, 25)
bimp = (0., 3.)


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    exe = pathlib.Path(EXE_DIR, path_main.stem)
    subprocess.run(
        f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(SRC_DIR)}', shell=True, text=True)

    path_data21 = pathlib.Path(
        input_dir, f'{reaction}E{energy}_{skyrme}_table21.root')
    path_data3 = pathlib.Path(
        input_dir, f'{reaction}E{energy}_{skyrme}_table3.root')

    path_out = pathlib.Path(
        out_dir, f'{reaction}E{energy}_{skyrme}_Ncmin{Nc[0]}_Ncmax{Nc[1]}.root')

    inputs = list(
        map(str, [f'{reaction}E{energy}', path_data21, path_data3, path_out, *Nc, *bimp]))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')


if __name__ == '__main__':
    main()
