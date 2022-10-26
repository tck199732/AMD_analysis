
import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyamd import PROJECT_DIR


class microball:
    def __init__(self, path=None, coordinate='uball'):
        if path is None:
            self.path = pathlib.Path(
                PROJECT_DIR, 'database/microball/acceptance/uball_geometry.dat')
        else:
            self.path = pathlib.Path(path)
        self.read_geometry()

        if coordinate == 'hira':
            self.hira_coordinate()

        self.read_e15190_config()

    def read_geometry(self):
        self.geometry = pd.read_csv(
            self.path, delim_whitespace=True, comment='#', index_col=False)
        self.geometry.set_index(['ring', 'det'], inplace=True)
        return self.geometry

    def hira_coordinate(self):

        self.geometry['phi_min'] = self.geometry['phi_min'] + 90.
        self.geometry['phi_max'] = self.geometry['phi_max'] + 90.

        for i, row in self.geometry.iterrows():
            print(i)
            phi_min = row['phi_min']
            phi_max = row['phi_max']
            if phi_min >= 360:
                self.geometry.at[i, 'phi_min'] -= 360
                #self.geometry.loc[i]['phi_min'] -= 360
            if phi_max >= 360:
                self.geometry.at[i, 'phi_max'] -= 360
            # elif phi_max < 0:
                #self.geometry.at[i, 'phi_max'] += 360

    def get_theta_range(self, ring, det):
        return self.geometry.loc[ring, det]['theta_min'], self.geometry.loc[ring, det]['theta_max']

    def get_phi_range(self, ring, det):
        return self.geometry.loc[ring, det]['phi_min'], self.geometry.loc[ring, det]['phi_max']

    def read_e15190_config(self, path=None, system=None):
        if path is None or system is None:
            self.config = {
                2: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                3: [1, 2,             7, 8, 9, 10, 11, 12],
                4: [1, 2, 3,          7, 8, 9, 10, 11, 12],
                5: [1, 2, 3,          7, 8, 9, 10, 11, 12, 13, 14],
                7: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                8: [1, 2, 3,    5, 6, 7, 8, 9, 10],
            }

        #path = "./acceptance/e15190_config.dat"
        # df = pd.read_csv(path, delim_whitespace = True, comment = '#', index_col = False)

    def draw_geometry(self, **kwargs):
        fig, ax = plt.subplots(**kwargs)
        cm = plt.cm.viridis(np.linspace(0.3, 1, len(self.config)))

        for ir, ring in enumerate(self.config):
            for det in self.config[ring]:
                theta_min, theta_max = self.get_theta_range(ring, det)
                phi_min, phi_max = self.get_phi_range(ring, det)
                if phi_max > phi_min:
                    box = plt.Rectangle(
                        (phi_min, theta_min),
                        phi_max - phi_min, theta_max - theta_min,
                        fill=True,
                        facecolor=cm[ir],
                        edgecolor='k',
                        linewidth=1,
                    )
                    ax.add_patch(box)

                    ax.text(
                        phi_min + 0.5 * (phi_max - phi_min),
                        theta_min + 0.5 * (theta_max - theta_min),
                        f'R{ring}-{det:02d}',
                        ha='center', va='center',
                    )
                else:
                    box_a = plt.Rectangle(
                        (phi_min, theta_min),
                        360. - phi_min, theta_max - theta_min,
                        fill=True,
                        facecolor=cm[ir],
                        edgecolor='k',
                        linewidth=1,
                    )
                    box_b = plt.Rectangle(
                        (0, theta_min),
                        phi_max, theta_max - theta_min,
                        fill=True,
                        facecolor=cm[ir],
                        edgecolor='k',
                        linewidth=1,
                    )
                    ax.add_patch(box_a)
                    ax.add_patch(box_b)

                    # ax.text(
                    #phi_min + 0.5 * (phi_max - phi_min),
                    #theta_min + 0.5 * (theta_max - theta_min),
                    # f'R{ring}-{det:02d}',
                    # ha='center', va='center',
                    # )

        ax.set_title('Microball coverage in E15190-E14030')
        ax.set_xlim(-20, 360)
        ax.set_ylim(0, 180)
        ax.set_xlabel(r'$\phi$ (deg)')
        ax.set_ylabel(r'$\theta$ (deg)')
        return fig, ax

    def impact_parameter_mapping(self, path=None):
        if path is None:
            path = pathlib.Path(
                PROJECT_DIR, 'database/microball/bimp_mapping/Ca48Ni64E140.dat')
        df = pd.read_csv(str(path), delim_whitespace=True)
        return df

    def plot_impact_parameter_mapping(self, ax=None, path=None, **kwargs):
        df = self.impact_parameter_mapping(path)
        if ax is None:
            ax = plt.gca()
        kw = dict(
            fmt='.'
        )
        kw.update(kwargs)
        ax.errorbar(df['multiplicity'], df['b'], yerr=df['b_err'], **kw)
        return ax


if __name__ == '__main__':

    uball = microball(
        "../../database/microball//acceptance/uball_geometry.dat", 'uball')
    # df = uball.geometry
    # print(df)
    # fig, ax = uball.draw_geometry(figsize=(10, 6))
    # fig.savefig('uball_coverage_e15190.png')
    fig, ax = plt.subplots()
    uball.plot_impact_parameter_mapping(ax)
    fig.savefig('impact_parameter_mapping.png')
