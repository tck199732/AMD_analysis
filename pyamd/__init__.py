import pathlib
import os
PROJECT_DIR = pathlib.Path(__file__).parent.parent.resolve()
DATA_DIR = pathlib.Path('/data/amd/dec2021/b3fm/filtered')
os.environ['PROJECT_DIR'] = str(PROJECT_DIR)
os.environ['DATA_DIR'] = str(DATA_DIR)
