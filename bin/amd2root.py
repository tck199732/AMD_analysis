import inspect
import subprocess
import pathlib
import os
import re
import glob

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, 'database')
dat_dir = pathlib.Path('/data/amd/dec2021/b3fm')
root_dir = pathlib.Path('/data/amd/dec2021/b3fm')
list_dir = pathlib.Path(database, 'inputlist/dec2021')

# inputs
beam = 'Ca'
target = 'Ni'
beamA = 48
targetA = 64
energy = 140
skyrme = 'SkM'

##########################################################################
reaction = f'{beam}{beamA}{target}{targetA}E{energy}'
rec_name = f'{beam}{beamA}{target}{targetA}En{energy}MeV_{skyrme}'
exe = pathlib.Path(str(project_dir), 'bin/amd2root')
path_list = pathlib.Path(list_dir, f'{rec_name}.list')


def main():
    while not exe.exists():
        try:
            os.chdir(f'{str(project_dir)}/bin')
            root_libs = subprocess.run('root-config --libs --glibs --cflags',
                                       shell=True, capture_output=True, text=True, encoding='utf-8')
            subprocess.run(
                f'g++ amd2root.cpp -o amd2root -I{root_libs.stdout.strip()}', shell=True)
        except:
            current_env = ''
            conda_envs = subprocess.run(
                'conda env list', capture_output=True, shell=True, text=True)
            envs = conda_envs.stdout.strip().split('\n')[2:]
            for env in envs:
                if '*' in env:
                    current_env = envs[0].split()[-1]

            if not current_env == f'{str(project_dir)}/env':
                os.chdir(str(project_dir))
                subprocess.run(f'conda activate ./env', shell=True)

    run(mode='21')
    # run(mode='21t')
    run(mode='3')
    print('All DONE')


def run(mode):

    m = re.sub('[a-z]', '', mode)
    with open(str(path_list), 'r') as file:
        fl = file.readlines()
        fl = [line for line in fl if not line.isspace()]
        path_data = [f'{str(dat_dir)}/{fn.strip()}_table{m}.dat' for fn in fl]
        path_out = [
            f'{str(root_dir)}/{fn.strip()}_table{mode}.root' for fn in fl]
        path_coll_hist = [
            f'{str(dat_dir)}/{fn.strip()}_coll_hist.dat' for fn in fl]
        path_amdgid = [f'{str(dat_dir)}/{fn.strip()}_amdgid.dat' for fn in fl]

    # print(path_data)
    # print(path_out)
    # print(path_coll_hist)
    # print(path_amdgid)

    for dat, out, ch, gid in zip(path_data, path_out, path_coll_hist, path_amdgid):
        inputs = [reaction, mode, dat, out]
        if mode == "21t":
            inputs.extend([ch, gid])
        args = ' '.join(inputs)
        subprocess.run(f'{str(exe)} {args}', shell=True, text=True)


if __name__ == '__main__':
    main()
