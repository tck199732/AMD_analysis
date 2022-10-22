import subprocess
import pathlib
import os
import itertools

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, f'database')
list_dir = pathlib.Path(database, 'inputlist/dec2021')
# list_dir = pathlib.Path(database, 'inputlist/feb2022')

exe_dir = pathlib.Path(project_dir, 'bin')
exe = pathlib.Path(exe_dir, 'ExpFilter')
src_dir = pathlib.Path(project_dir, 'src')
input_dir = pathlib.Path('/data/amd/dec2021/b3fm')
# input_dir = pathlib.Path('/data/amd/feb2022/b10fm')
out_dir = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
# out_dir = pathlib.Path(f'/data/amd/dec2022/b10fm/filtered')
out_dir.mkdir(exist_ok=True)


nuclei = ['Ca40Ni58', 'Ca48Ni64', 'Ca40Sn112', 'Ca48Sn124']
energy = [56, 140]
skyrme = ['SkM', 'SLy4', 'SLy4_L108']
combination = list(itertools.product(nuclei, energy, skyrme))
reaction = [(f'{nuc}E{e}', sky) for nuc, e, sky in combination]
rec_name = [f'{nuc}En{e}MeV_{sky}' for nuc, e, sky in combination]
path_list = {rec: pathlib.Path(list_dir, f'{name}.list')
             for rec, name in zip(reaction, rec_name)}


def main():

    if not exe.exists():
        os.chdir(exe_dir)
        subprocess.run(
            f'g++ ExpFilter.cpp -o ExpFilter -I`root-config --cflags --libs --glibs` -I{str(src_dir)}', shell=True)

    for i, (rec, rlist) in enumerate(path_list.items()):
        for mode in ['21', '3']:
            path_out = pathlib.Path(out_dir, f'{rec[0]}_{rec[1]}.root')
            inputs = list(map(str, [rec[0], mode, input_dir, rlist, path_out]))
            args = ' '.join(inputs)
            # print(f'{str(exe)} {args}')
            subprocess.run(f'{str(exe)} {args}', shell=True, text=True)

    print('All DONE')


if __name__ == '__main__':
    main()
