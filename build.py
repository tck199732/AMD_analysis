#!/usr/bin/env python
import sys
if not (sys.version_info[0] >= 3 and sys.version_info[1] >= 7):
    raise Exception('Requires at least Python 3.7 to use this repository.')

import pathlib
import subprocess
import inspect
import glob

# env_name is used only as the directory name of the env 
env_name = 'env_amd' 
project_dir = pathlib.Path(__file__).parent.resolve()
env_yml = pathlib.Path(project_dir, 'environment.yml')
env_dir = pathlib.Path(project_dir, env_name)

def main():
    
    check_conda_activated()

    # create conda environment
    subprocess.run(f'conda env create -f {env_yml} --prefix {env_dir}', shell=True, text=True)
    
    # create script to save environment variables
    env_vars = {
        'PROJECT_DIR' : str(project_dir),
    }

    scripts = {

    'activate' : [inspect.cleandoc(f"""
        if [ ! -z '${name}' ]; then
            export OLD_{name}=${name};
        fi
        export {name}='{value}'
    """) for name, value in env_vars.items()],
    
    'deactivate' : [inspect.cleandoc(f"""
        if [ ! -z '$OLD_{name}' ]; then
            {name}='$OLD_{name}'
            unset OLD_{name}
        fi
    """) for name, value in env_vars.items()]

    }
    
    local_dir = pathlib.Path(project_dir, 'local')
    for name in ['activate', 'deactivate']:

        src_dir = pathlib.Path(local_dir, f'etc/conda/{name}.d')
        src_dir.mkdir(exist_ok=True, parents=True)
        content = '\n'.join(['#!/bin/bash'] + scripts[name])
        src_pth = pathlib.Path(src_dir, 'env_vars.sh')

        with open(str(src_pth), 'w') as f:
            f.write(content + '\n')
        
        target_dir =  pathlib.Path(env_dir, f'etc/conda/{name}.d')
        target_dir.mkdir(exist_ok=True, parents=True)
        target_pth = pathlib.Path(target_dir, 'env_vars.sh')
        if target_pth.is_symlink():
            target_pth.unlink()
        target_pth.symlink_to(src_pth)
        
    # add the project home directory to conda.pth
    pths = glob.glob(f'{str(env_dir)}/lib/python*/site-packages')
    for pth in pths:
        pth = pathlib.Path(pth, 'conda.pth')
        with open(str(pth), 'w') as f:
            f.write(str(project_dir) + '\n')


def check_conda_activated():
    which_conda = subprocess.run('which conda', shell=True, capture_output=True, text=True)
    path1 = pathlib.Path(which_conda.stdout.strip())
    path1 = pathlib.Path(path1.parent.parent, 'conda-meta')
    path2 = pathlib.Path(sys.prefix, 'conda-meta')
    if not (path1.is_dir() and path2.is_dir()):
        raise Exception('Fail to detect conda environment. Please make sure you have activated conda.')

if __name__ == '__main__':
    main()