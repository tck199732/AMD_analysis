import pathlib
import os
PROJECT_DIR = pathlib.Path(__file__).parent.parent.resolve()
os.environ['PROJECT_DIR'] = str(PROJECT_DIR)
