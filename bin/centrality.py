from pyamd import PROJECT_DIR, DATABASE, EXE_DIR, SRC_DIR
import subprocess
import pathlib
import os
import time


exe_name = pathlib.Path(__file__).stem
path_main = pathlib.Path(EXE_DIR, f'{str(exe_name)}.cpp')
exe = pathlib.Path(EXE_DIR, exe_name)

# input_dir = pathlib.Path('/data/amd/feb2022/b10fm/filtered')
# input_dir = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
input_dir = pathlib.Path('/data/amd/nov2022/sigma0.85/filtered')

# out_dir = pathlib.Path(PROJECT_DIR, f'result/{exe_name}/b10fm')
out_dir = pathlib.Path(PROJECT_DIR, f'result/{exe_name}/b3fm')
out_dir.mkdir(exist_ok=True, parents=True)

reaction = 'Ca48Sn124'
energy = 140
mode = '3'
skyrme = 'SkM'


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    subprocess.run(
        f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(SRC_DIR)}', shell=True, text=True)

    path_data = pathlib.Path(
        input_dir, f'{reaction}E{energy}_{skyrme}_table{mode}.root')
    path_out = pathlib.Path(
        out_dir, f'{reaction}E{energy}_{skyrme}_table{mode}.root')
    inputs = list(
        map(str, [f'{reaction}E{energy}', mode, path_data, path_out]))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'time elapsed = {elapsed_time:.2f}')


if __name__ == '__main__':
    main()
