import pathlib
import os
PROJECT_DIR = pathlib.Path(os.environ['CONDA_PREFIX']).parent
os.environ['PROJECT_DIR'] = str(PROJECT_DIR)
DATABASE = pathlib.Path(PROJECT_DIR, f'database')
EXE_DIR = pathlib.Path(PROJECT_DIR, 'bin')
ANAL_DIR = pathlib.Path(PROJECT_DIR, 'analysis')
SRC_DIR = pathlib.Path(PROJECT_DIR, 'src')
