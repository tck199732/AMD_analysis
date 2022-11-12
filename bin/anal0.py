import subprocess
import pathlib
import os
import itertools
import time

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, f'database')

#######################################################################
exe_dir = pathlib.Path(project_dir, 'bin')
exe = pathlib.Path(exe_dir, 'anal0')
src_dir = pathlib.Path(project_dir, 'src')

list_dir = pathlib.Path(database, 'inputlist/dec2021')
input_dir = pathlib.Path('/data/amd/dec2021/b3fm')
out_dir = pathlib.Path(project_dir, 'result/anal0/sigma0.85_b3fm')

# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma100')
# input_dir = pathlib.Path('/data/amd/nov2022/sigma100')
# out_dir = pathlib.Path(project_dir, 'result/anal0/sigma100_b3fm')

# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma_free')
# input_dir = pathlib.Path('/data/amd/nov2022/sigma_free')
# out_dir = pathlib.Path(project_dir, 'result/anal0/sigma_free_b3fm')
#######################################################################
out_dir.mkdir(exist_ok=True, parents=True)



REACTION = 'Ca40Ni58E56'
SKYRME = 'SkM'

# a collection of all reaction systems, only run program if path exists
nuclei = ['Ca40Ni58', 'Ca48Ni64', 'Ca40Sn112', 'Ca48Sn124']
energy = [56, 140]
skyrme = ['SkM', 'SLy4', 'SLy4_L108']
combination = list(itertools.product(nuclei, energy, skyrme))
reaction = [(f'{nuc}E{e}', sky) for nuc, e, sky in combination]
rec_name = [f'{nuc}En{e}MeV_{sky}' for nuc, e, sky in combination]
path_list = {rec: pathlib.Path(list_dir, f'{name}.list')
             for rec, name in zip(reaction, rec_name)}


def main():
    os.chdir(exe_dir)
    subprocess.run(
        f'g++ anal0.cpp -o anal0 -I`root-config --cflags --libs --glibs` -I{str(src_dir)}', shell=True)

    run(REACTION, SKYRME)
    print('All DONE')


def run(reaction, skyrme):
    pth_ls = path_list[(reaction, skyrme)]
    if not pth_ls.exists():
        raise OSError(f'{str(pth_ls)} does not exist.')
    path_out = pathlib.Path(
        out_dir, f'{reaction}_{skyrme}.root')
    inputs = list(map(str, [reaction, input_dir, pth_ls, path_out]))
    args = ' '.join(inputs)
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("elapsed time = ", elapsed_time)


if __name__ == '__main__':
    main()
