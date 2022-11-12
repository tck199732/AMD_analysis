import subprocess
import pathlib
import os
import time

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, f'database')

exe_dir = pathlib.Path(project_dir, 'bin')
src_dir = pathlib.Path(project_dir, 'src')
exe_name = pathlib.Path(__file__).stem
path_main = pathlib.Path(exe_dir, f'{str(exe_name)}.cpp')

# input_dir = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
# out_dir = pathlib.Path(project_dir, f'result/{str(exe_name)}/sigma0.85_b3fm')

# input_dir = pathlib.Path('/data/amd/nov2022/sigma100/filtered')
# out_dir = pathlib.Path(project_dir, f'result/{str(exe_name)}/sigma100_b3fm')

input_dir = pathlib.Path('/data/amd/nov2022/sigma_free/filtered')
out_dir = pathlib.Path(project_dir, f'result/{str(exe_name)}/sigma_free_b3fm')

out_dir.mkdir(exist_ok=True, parents=True)

reaction = 'Ca40Ni58'
energy = 56
skyrme = 'SkM'
mode = "3"
Nc = (11, 25)
bimp = (0., 3.)


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    exe = pathlib.Path(exe_dir, path_main.stem)
    subprocess.run(
        f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(src_dir)}', shell=True, text=True)

    path_data = pathlib.Path(
        input_dir, f'{reaction}E{energy}_{skyrme}_table{mode}.root')

    path_out = pathlib.Path(
        out_dir, f'{reaction}E{energy}_{skyrme}_Ncmin{Nc[0]}_Ncmax{Nc[1]}_table{mode}.root')

    inputs = list(
        map(str, [f'{reaction}E{energy}', mode, path_data, path_out, *Nc, *bimp]))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')


if __name__ == '__main__':
    main()
