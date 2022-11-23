from pyamd import PROJECT_DIR, SRC_DIR, EXE_DIR, DATABASE
import subprocess
import pathlib
import os
import time

# project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
# database = pathlib.Path(project_dir, f'database')
# exe_dir = pathlib.Path(project_dir, 'bin')
# src_dir = pathlib.Path(project_dir, 'src')

path_main = pathlib.Path(EXE_DIR, 'spacetime.cpp')
input_dir = pathlib.Path('/data/amd/nov2022/sigma100')
out_dir = pathlib.Path(PROJECT_DIR, 'result/spacetime/sigma100_b3fm')

out_dir.mkdir(exist_ok=True, parents=True)
reaction = 'Ca40Ni58'
energy = 140
skyrme = 'SkM'


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    exe = pathlib.Path(EXE_DIR, path_main.stem)
    subprocess.run(
        f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(SRC_DIR)}', shell=True, text=True)

    path_data21t = pathlib.Path(
        input_dir, f'{reaction}_En{energy}MeV_{skyrme}_33k_table21t.root')

    path_out = pathlib.Path(out_dir, f'{reaction}E{energy}_{skyrme}.root')

    inputs = list(
        map(str, [f'{reaction}E{energy}', path_data21t, path_out]))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')


if __name__ == '__main__':
    main()
