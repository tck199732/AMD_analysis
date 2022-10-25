
import re
import numpy as np


class particle:
    PARTICLES = ['n', 'p', 'd', 't', '3He', '4He']
    ZID = [0, 1, 1, 1, 2, 2]
    NID = [1, 0, 1, 2, 1, 2]

    def __init__(self, name):
        if name == 'coal_p':
            name = 'p'
        elif name == 'coal_n':
            name = 'n'
        elif not name in self.PARTICLES:
            raise ValueError('input particle is not analyzed in experiment.')
        self.name = name
        self.Z = self.ZID[self.PARTICLES.index(name)]
        self.N = self.NID[self.PARTICLES.index(name)]


class reaction:

    mass_1u = 931.49410242
    beam_mass = {
        'Ca40': 39.962590866,
        'Ca48': 47.95252290
    }
    target_mass = {
        'Ni58': 57.935343,
        'Ni64': 63.9279660,
        'Sn112': 111.904818,
        'Sn124': 123.9052739
    }

    def __init__(self, system):
        self.beamA, self.targetA, self.beam_energy = [
            int(m) for m in re.compile('[0-9]+').findall(system)]
        self.beam, self.target = [
            m for m in re.compile('[A-Za-z]{2}').findall(system)]

        self.beam_mass = (
            self.beam_mass[f'{self.beam}{self.beamA}']) * self.mass_1u
        self.target_mass = (
            self.target_mass[f'{self.target}{self.targetA}']) * self.mass_1u

    def get_betacms(self):
        beam_ke = self.beam_energy * self.beamA
        beam_energy_tot = beam_ke + self.beam_mass
        mom_beam = np.sqrt(beam_ke**2 + 2*beam_ke * self.beam_mass)

        gamma = beam_energy_tot/self.beam_mass
        return mom_beam/(gamma*self.beam_mass + self.target_mass)

    def get_rapidity_beam(self):

        beam_ke = self.beam_energy * self.beamA
        mom_beam = np.sqrt(beam_ke**2 + 2*beam_ke * self.beam_mass)
        return 0.5 * np.log((beam_ke + self.beam_mass + mom_beam) / (beam_ke + self.beam_mass - mom_beam))
        # return get_rapidity(beam_ke, self.beam_mass, mom_beam)
