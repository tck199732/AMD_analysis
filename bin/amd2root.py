import inspect
import subprocess
import pathlib
import os
import re
import glob

home_dir = pathlib.Path(os.environ['HOME'])
project_dir = pathlib.Path(__file__).parent.parent.resolve()

dat_dir = pathlib.Path('/data/amd/dec2021')
root_dir = pathlib.Path('/data/amd/dec2021')

# inputs
beam = 'Ca'
target = 'Ni'
beamA = 48
targetA = 64
energy = 140
skyrme = 'SkM'

##########################################################################
reaction = f'{beam}{beamA}{target}{targetA}E{energy}'
exe = pathlib.Path(str(project_dir), 'bin/amd2root')


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
    run(mode='21t')
    run(mode='3')
    print('All DONE')


def run(mode):

    m = re.sub('[a-z]', '', mode)
    path_data = glob.glob(f'{str(dat_dir)}/table{m}.dat')
    path_out = glob.glob(f'{str(dat_dir)}/table{mode}.root')
    path_coll_hist = glob.glob(f'{str(dat_dir)}/coll_hist.dat')
    path_amdgid = glob.glob(f'{str(dat_dir)}/amdgid.dat')

    path_data = sorted([f for f in path_data if reaction.lower() in f.lower()])
    path_out = sorted([f for f in path_out if reaction.lower() in f.lower()])
    path_coll_hist = sorted(
        [f for f in path_coll_hist if reaction.lower() in f.lower()])
    path_amdgid = sorted(
        [f for f in path_amdgid if reaction.lower() in f.lower()])

    os.chdir(f'{str(project_dir)}/bin')

    for dat, out, ch, gid in zip(path_data, path_out, path_coll_hist, path_amdgid):
        inputs = [reaction, mode, dat, out]
        if mode == "21t":
            inputs.extend([ch, gid])
        args = ' '.join(inputs)
        subprocess.run(f'{str(exe)} {args}', shell=True, text=True)


if __name__ == '__main__':
    main()
