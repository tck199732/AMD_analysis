import subprocess
import pathlib
import os


project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, f'database')

src_dir = pathlib.Path(project_dir, 'src')
exe_dir = pathlib.Path(project_dir, 'bin')
path_main = pathlib.Path(exe_dir, 'multiplicity.cpp')
exe = pathlib.Path(exe_dir, path_main.stem)

input_dir = pathlib.Path('/data/amd/feb2022/b10fm/filtered')
out_dir = pathlib.Path(project_dir, 'result/multiplicity')
out_dir.mkdir(exist_ok=True)

reaction = 'Ca48Ni64'
energy = 140
mode = '3'
skyrme = 'SkM'


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs', capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()

    subprocess.run(f'g++ {str(path_main)} -o {str(exe)} -I{root_inc} -I{str(src_dir)}', shell=True, text=True)

    path_data = pathlib.Path(
        input_dir, f'{reaction}E{energy}_{skyrme}_table{mode}.root')
    path_out = pathlib.Path(
        out_dir, f'{reaction}E{energy}_{skyrme}_table{mode}_b10fm.root')
    inputs = list(
        map(str, [f'{reaction}E{energy}', mode, path_data, path_out]))
    args = ' '.join(inputs)
    # print(f'{str(exe)} {args}')
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)


if __name__ == '__main__':
    main()
