from pickle import TRUE
import subprocess
import pathlib
import os

env_dir = pathlib.Path(os.environ['CONDA_PREFIX'])
project_dir = env_dir.parent.resolve()
exe_dir = pathlib.Path(project_dir, 'bin')
# out_dir = pathlib.Path(project_dir, 'result')
out_dir = pathlib.Path('/data/AMD/filtered')
out_dir.mkdir(exist_ok=True)

beam = 'Ca'
beamA = 40
target = 'Ni'
targetA = 58
energy = 140
skyrme = 'SkM'


def main():
    reaction = f'{beam}{beamA}{target}{targetA}E{energy}'
    mode = '21'
    path_list = pathlib.Path(
        exe_dir, f'filelist/dec2021/{beam}{beamA}{target}{targetA}_En{energy}MeV_{skyrme}.list')

    path_out = pathlib.Path(
        out_dir, f'{beam}{beamA}{target}{targetA}_En{energy}MeV_{skyrme}_table{mode}.root')

    os.chdir(exe_dir)

    subprocess.run(
        f'a {reaction} {mode} {path_list} {path_out}', shell=True, text=True)


if __name__ == '__main__':
    main()
