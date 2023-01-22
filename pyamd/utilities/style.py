import json
import pathlib
from pyamd import DATABASE

import pathlib
import json

PATH_STYLE = pathlib.Path(DATABASE, 'styles/matplotlib.json')


def export_mpl_style(path=PATH_STYLE, mode='w'):
    dir = pathlib.Path(path).parent.resolve()
    dir.mkdir(exist_ok=True, parents=True)

    style = {
        "figure.facecolor": "white",
        "figure.edgecolor": "black",
        "figure.frameon": True,
        "figure.dpi": 150,
        "figure.figsize": [4, 3.5],

        "axes.facecolor":     "white",   # axes background color
        "axes.edgecolor":     "black",  # axes edge color
        "axes.linewidth":     0.8,     # edge line width
        "axes.grid":          False,   # display grid or not
        "axes.grid.axis":     "both",    # which axis the grid should apply to
        # grid lines at {major, minor, both} ticks
        "axes.grid.which":    "major",
        # alignment of the title: {left, right, center}
        "axes.titlelocation": "center",
        "axes.titlesize":     15,   # font size of the axes title
        "axes.titleweight":   "normal",  # font weight of title
        "axes.titlecolor":    "auto",    # color of the axes title, auto falls back to
        # text.color as default value
        # position title (axes relative units).  None implies auto
        "axes.titley":        "None",
        "axes.titlepad":      6.0,     # pad between axes and title in points
        "axes.labelsize":     12.5,  # font size of the x and y labels
        "axes.labelpad":      4.0,     # space between label and axis
        "axes.labelweight":   "normal",  # weight of the x and y labels
        "axes.labelcolor":    "black",
        "axes.axisbelow":     "line",    # draw axis gridlines and ticks:
        #     - below patches (True)
        #     - above patches but below lines ('line')
        #     - above all (False)

        "font.family": "serif",
        "font.size": 15,

        "mathtext.fontset": "cm",
        # "mathtext.bf" :,
        # "mathtext.cal" :,
        # "mathtext.default" :,
        # "mathtext.fallback" :,
        # "mathtext.fontset" :,
        # "mathtext.it" :,
        # "mathtext.rm" :,
        # "mathtext.sf" :,
        # "mathtext.tt" :,

        "lines.color": "black",
        "lines.linestyle": '--',
        "lines.linewidth": 2.,

        "grid.linestyle": "dashed",
        "grid.alpha": 0.8,
        "grid.color": "black",
        "grid.linewidth": 1.,

        # ***************************************************************************
        # * TICKS                                                                   *
        # ***************************************************************************

        "xtick.top":           True,  # draw ticks on the top side
        "xtick.bottom":        True,  # draw ticks on the bottom side
        "xtick.labeltop":      False,  # draw label on the top
        "xtick.labelbottom":   True,  # draw label on the bottom
        "xtick.major.size":    3.5,  # major tick size in points
        "xtick.minor.size":    2,  # minor tick size in points
        "xtick.major.width":   0.8,  # major tick width in points
        "xtick.minor.width":   0.6,  # minor tick width in points
        "xtick.major.pad":     3.5,  # distance to major tick label in points
        "xtick.minor.pad":     3.4,  # distance to the minor tick label in points
        "xtick.color":         "black",  # color of the ticks
        # color of the tick labels or inherit from xtick.color
        "xtick.labelcolor":    "inherit",
        "xtick.labelsize":     "small",  # font size of the tick labels
        "xtick.direction":     "in",  # direction: {in, out, inout}
        "xtick.minor.visible": True,  # visibility of minor ticks on x-axis
        "xtick.major.top":     True,  # draw x axis top major ticks
        "xtick.major.bottom":  True,  # draw x axis bottom major ticks
        "xtick.minor.top":     True,  # draw x axis top minor ticks
        "xtick.minor.bottom":  True,  # draw x axis bottom minor ticks
        "xtick.alignment":     "center",  # alignment of xticks

        "ytick.left":          True,  # draw ticks on the left side
        "ytick.right":         True,  # draw ticks on the right side
        "ytick.labelleft":     True,  # draw tick labels on the left side
        "ytick.labelright":    False,  # draw tick labels on the right side
        "ytick.major.size":    3.5,  # major tick size in points
        "ytick.minor.size":    2,  # minor tick size in points
        "ytick.major.width":   0.8,  # major tick width in points
        "ytick.minor.width":   0.6,  # minor tick width in points
        "ytick.major.pad":     3.5,  # distance to major tick label in points
        "ytick.minor.pad":     3.4,  # distance to the minor tick label in points
        "ytick.color":         "black",  # color of the ticks
        # color of the tick labels or inherit from ytick.color
        "ytick.labelcolor":    "inherit",
        "ytick.labelsize":     "small",  # font size of the tick labels
        "ytick.direction":     "in",  # direction: {in, out, inout}
        "ytick.minor.visible": True,  # visibility of minor ticks on y-axis
        "ytick.major.left":    True,  # draw y axis left major ticks
        "ytick.major.right":   True,  # draw y axis right major ticks
        "ytick.minor.left":    True,  # draw y axis left minor ticks
        "ytick.minor.right":   True,  # draw y axis right minor ticks
        "ytick.alignment":     "center_baseline",   # alignment of yticks

        # ***************************************************************************
        # * LEGEND                                                                   *
        # ***************************************************************************

        "legend.loc":           "best",
        "legend.frameon":       True,     # if True, draw the legend on a background patch
        "legend.framealpha":    0.8,      # legend patch transparency
        "legend.facecolor":     "inherit",  # inherit from axes.facecolor; or color spec
        "legend.edgecolor":     "None",      # background patch boundary color
        # if True, use a rounded box for the legend background, else a rectangle
        "legend.fancybox":      True,
        "legend.shadow":        False,    # if True, give background a shadow effect
        "legend.numpoints":     1,        # the number of marker points in the legend line
        "legend.scatterpoints": 1,        # number of scatter points
        "legend.markerscale":   1.0,      # the relative size of legend markers vs. original
        "legend.fontsize":      "x-small",
        "legend.labelcolor":    "None",
        # None sets to the same as the default axes.
        "legend.title_fontsize": "None",

        # Dimensions as fraction of font size:
        "legend.borderpad":     0.3,  # border whitespace
        "legend.labelspacing":  0.5,  # the vertical space between the legend entries
        "legend.handlelength":  2.0,  # the length of the legend lines
        "legend.handleheight":  0.7,  # the height of the legend handle
        "legend.handletextpad": 0.8,  # the space between the legend line and legend text
        "legend.borderaxespad": 0.5,  # the border between the axes and legend edge
        "legend.columnspacing": 2.0,  # column separation
    }
    with open(str(path), mode) as f:
        json.dump(style, f)


def set_matplotlib_style(mpl):
    with open(PATH_STYLE, 'r') as f:
        style = json.load(f)
    mpl.rcParams.update(style)


if __name__ == '__main__':
    export_mpl_style()