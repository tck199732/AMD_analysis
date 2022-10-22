import subprocess
import pathlib
import os
import itertools

anal_date = 'dec2021'
env_dir = pathlib.Path(os.environ['CONDA_PREFIX'])
project_dir = env_dir.parent.resolve()
database = pathlib.Path(project_dir, f'database')
exe_dir = pathlib.Path(project_dir, 'bin')
src_dir = pathlib.Path(project_dir, 'src')
out_dir = pathlib.Path(f'/data/AMD/filtered/{anal_date}')
out_dir.mkdir(exist_ok=True)

beam = 'Ca'
beamA = 48
target = 'Ni'
targetA = 64
energy = 140
skyrme = 'SkM'
mode = "21"
mode = "3"
##############################################################################
# beam = ['Ca40', 'Ca48']
# target = ['Ni58', 'Ni64', 'Sn112', 'Sn124']
# energy = [56,140]
# mode = ["21", "3"]
##############################################################################


def main():
    reaction = f'{beam}{beamA}{target}{targetA}E{energy}'
    path_list = pathlib.Path(
        exe_dir, f'{str(database)}/inputlist/dec2021/{beam}{beamA}{target}{targetA}En{energy}MeV_{skyrme}.list')

    path_out = pathlib.Path(
        out_dir, f'{beam}{beamA}{target}{targetA}_En{energy}MeV_{skyrme}_table{mode}.root')

    filter = pathlib.Path(exe_dir, 'ExpFilter')
    if not filter.exists():
        os.chdir(exe_dir)
        subprocess.run(
            f'g++ ExpFilter.cpp -o ExpFilter -I`root-config --cflags --libs --glibs` -I{str(src_dir)}', shell=True)

    print('reaction : ', reaction)
    print('mode : ', mode)
    print(str(path_list))
    print(str(path_out))
    print(str(filter))

    subprocess.run(
        f'{str(filter)} {reaction} {mode} {path_list} {path_out}', shell=True, text=True)


if __name__ == '__main__':
    main()
