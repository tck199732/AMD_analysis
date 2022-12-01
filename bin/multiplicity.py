from pyamd import PROJECT_DIR, DATABASE, EXE_DIR, SRC_DIR
import subprocess
import pathlib
import os
import time


exe_name = pathlib.Path(__file__).stem
path_main = pathlib.Path(EXE_DIR, f'{str(exe_name)}.cpp')

input_dir = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
out_dir = pathlib.Path(PROJECT_DIR, f'result/{str(exe_name)}/sigma0.85_b3fm')
# input_dir = pathlib.Path('/data/amd/nov2022/sigma100/filtered')
# out_dir = pathlib.Path(PROJECT_DIR, f'result/{str(exe_name)}/sigma100_b3fm')
out_dir.mkdir(exist_ok=True, parents=True)

reaction = 'Ca40Ni58'
energy = 140
skyrme = 'SkM'
mode = "3"
Nc = (11, 25)
bimp = (0., 3.)


def main():

    root_inc = 'root-config --cflags --libs --glibs'
    root_inc = subprocess.run(root_inc, capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    exe = pathlib.Path(EXE_DIR, path_main.stem)
    exe_cmd = f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(SRC_DIR)}'
    subprocess.run(exe_cmd, shell=True, text=True)

    data_fname = f'{reaction}E{energy}_{skyrme}_table{mode}.root'
    path_data = pathlib.Path(input_dir, data_fname)

    out_fname = f'{reaction}E{energy}_{skyrme}_Ncmin{Nc[0]}_Ncmax{Nc[1]}_table{mode}.root'
    path_out = pathlib.Path(out_dir, out_fname)


    inputs = [f'{reaction}E{energy}', mode, path_data, path_out, *Nc, *bimp]
    inputs = list(map(str, inputs))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')

if __name__ == '__main__':
    main()
