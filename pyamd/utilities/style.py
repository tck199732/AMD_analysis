import json
import pathlib

from pyamd import PROJECT_DIR
def set_matplotlib_style(mpl):
    path = pathlib.Path(PROJECT_DIR, 'database/styles/matplotlib.json')
    with open(path, 'r') as file:
        content = json.load(file)
    mpl.rcParams.update(content)
        
