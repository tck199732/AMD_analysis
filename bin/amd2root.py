import inspect
import itertools
import subprocess
import pathlib
import os
import re
import glob

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, 'database')

# dat_dir = pathlib.Path('/data/amd/dec2021/b3fm')
# root_dir = pathlib.Path('/data/amd/dec2021/b3fm')
# list_dir = pathlib.Path(database, 'inputlist/dec2021')

# dat_dir = pathlib.Path('/data/amd/feb2022/b10fm')
# root_dir = pathlib.Path('/data/amd/feb2022/b10fm')
# list_dir = pathlib.Path(database, 'inputlist/feb2022')

# dat_dir = pathlib.Path('/data/amd/nov2022/sigma100')
# root_dir = pathlib.Path('/data/amd/nov2022/sigma100')
# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma100')

# dat_dir = pathlib.Path('/data/amd/nov2022/sigma_free')
# root_dir = pathlib.Path('/data/amd/nov2022/sigma_free')
# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma_free')

# dat_dir = pathlib.Path('/data/amd/nov2022/sigma10.05')
# root_dir = pathlib.Path('/data/amd/nov2022/sigma10.05')
# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma10.05')

dat_dir = pathlib.Path('/data/amd/nov2022/sigma0.25')
root_dir = pathlib.Path('/data/amd/nov2022/sigma0.25')
list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma0.25')

# dat_dir = pathlib.Path('/data/amd/nov2022/sigma5.05')
# root_dir = pathlib.Path('/data/amd/nov2022/sigma5.05')
# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma5.05')

# dat_dir = pathlib.Path('/data/amd/nov2022/sigma0.45')
# root_dir = pathlib.Path('/data/amd/nov2022/sigma0.45')
# list_dir = pathlib.Path(database, 'inputlist/nov2022/sigma0.45')

##########################################################################
nuclei = ['Ca40Ni58', 'Ca48Ni64', 'Ca40Sn112', 'Ca48Sn124']
energy = [56, 140]
skyrme = ['SkM', 'SLy4', 'SLy4_L108']
combination = list(itertools.product(nuclei, energy, skyrme))
reaction = [(f'{nuc}E{e}', sky) for nuc, e, sky in combination]
rec_name = [f'{nuc}En{e}MeV_{sky}' for nuc, e, sky in combination]
##########################################################################
exe = pathlib.Path(str(project_dir), 'bin/amd2root')
path_list = {rec: pathlib.Path(list_dir, f'{name}.list')
             for rec, name in zip(reaction, rec_name)}


def main():
    os.chdir(f'{str(project_dir)}/bin')
    root_libs = subprocess.run('root-config --libs --glibs --cflags',
                               shell=True, capture_output=True, text=True, encoding='utf-8')
    subprocess.run(
        f'g++ amd2root.cpp -o amd2root -I{root_libs.stdout.strip()}', shell=True)

    # while not exe.exists():
    #     try:
    #         os.chdir(f'{str(project_dir)}/bin')
    #         root_libs = subprocess.run('root-config --libs --glibs --cflags',
    #                                    shell=True, capture_output=True, text=True, encoding='utf-8')
    #         subprocess.run(
    #             f'g++ amd2root.cpp -o amd2root -I{root_libs.stdout.strip()}', shell=True)
    #     except:
    #         current_env = ''
    #         conda_envs = subprocess.run(
    #             'conda env list', capture_output=True, shell=True, text=True)
    #         envs = conda_envs.stdout.strip().split('\n')[2:]
    #         for env in envs:
    #             if '*' in env:
    #                 current_env = envs[0].split()[-1]

    #         if not current_env == f'{str(project_dir)}/env':
    #             os.chdir(str(project_dir))
    #             subprocess.run(f'conda activate ./env', shell=True)

    run(mode='21')
    run(mode='21t')
    run(mode='3')
    print('All DONE')


def run(mode):

    m = re.sub('[a-z]', '', mode)
    for i, (rec, flist) in enumerate(path_list.items()):
        if not flist.exists():
            continue
        with open(str(flist), 'r') as file:
            fl = file.readlines()
            fl = [line for line in fl if not line.isspace()]
            path_data = [
                pathlib.Path(f'{str(dat_dir)}/{fn.strip()}_table{m}.dat') for fn in fl]
            path_out = [
                pathlib.Path(f'{str(root_dir)}/{fn.strip()}_table{mode}.root') for fn in fl]
            path_coll_hist = [
                pathlib.Path(f'{str(dat_dir)}/{fn.strip()}_coll_hist.dat') for fn in fl]
            path_amdgid = [
                pathlib.Path(f'{str(dat_dir)}/{fn.strip()}_amdgid.dat') for fn in fl]

            # print(list(map(str, path_data)))
            # print(list(map(str, path_out)))
            # print(list(map(str, path_coll_hist)))
            # print(list(map(str, path_amdgid)))

        for dat, out, ch, gid in zip(path_data, path_out, path_amdgid, path_coll_hist):
            if out.exists():
                continue

            inputs = [rec[0], mode, dat, out]
            inputs = list(map(str, inputs))
            if mode == '21t':
                inputs.extend(list(map(str, [ch, gid])))
            args = ' '.join(inputs)
            print(args)
            subprocess.run(f'{str(exe)} {args}', shell=True, text=True)


if __name__ == '__main__':
    main()
