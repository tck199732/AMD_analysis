import pathlib
import re
import numpy as np
import itertools
from astropy import units
from collections import defaultdict

from pyamd import PROJECT_DIR
from pyamd.utilities import ame
ame_table = ame.AME()

class RunLog:
    CONFIG = pathlib.Path(PROJECT_DIR, 'database/e15190/RunInfo.data')
    def __init__(self, path_config=None):
        if path_config is None:
            path_config = self.CONFIG
        
        self.path_config = path_config

        self.RunIndices = dict()
        self.TriggerConditions = defaultdict(dict)
        self.BadMapVersion = defaultdict(dict)
        self.ShadowBar = defaultdict(dict)
        self._read()
        

    def _read(self):
        if not self.path_config.exists():
            raise OSError(f'{str(self.path_config)} does not exist.')

        with open(str(self.path_config), 'r') as f:
            content = f.readlines()
            for batch in content:
                if batch.startswith('#') or batch.isspace():
                    continue

                info = batch.split()
                reaction = info[0]
                startID, endID = int(info[1]), int(info[2])
                badmap = info[3]
                shadowbar = int(info[4])
                trigger = info[5]

                if not reaction in self.RunIndices:
                    self.RunIndices[reaction] = []
                
                for id in range(startID, endID+1):
                    self.RunIndices[reaction].append(id)
                    self.TriggerConditions[reaction][id] = trigger
                    self.BadMapVersion[reaction][id] = badmap
                    self.ShadowBar[reaction][id] = shadowbar


class Reaction:
    '''
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
    '''
    BEAM = ['Ca40', 'Ca48']
    TARGET = ['Ni58', 'Ni64', 'Sn112', 'Sn124']
    ENERGY = [56, 140]

    def __init__(self):
        self.reactions = [f'{b}{t}E{e}'.lower() for b, t, e in itertools.product(self.BEAM, self.TARGET, self.ENERGY)]
        self._initialize_reaction = False

    def set_reaction(self, reaction):
        if not reaction.lower() in self.reactions:
            raise ValueError('Invalid input name for reaction.')
        
        self._initialize_reaction = True

        self.beamA, self.targetA, beam_energy = [
            int(m) for m in re.compile('[0-9]+').findall(reaction)]
        self.beam, self.target = [
            m for m in re.compile('[A-Za-z]{2}').findall(reaction)]

        self.beam_energy = beam_energy * units.MeV
        beam = f'{self.beam}{self.beamA}'
        target = f'{self.target}{self.targetA}'

        self.beam_mass = ame_table.get_mass(beam)
        self.target_mass = ame_table.get_mass(target)
        

    def get_name(self, latex=False):
        if not self._initialize_reaction:
            raise ValueError('Choose the reaction before parsing its symbol.')

        if latex:
            name = r'$^{{beamA}}{beam} + ^{{targetA}}{target} @ {beamE}$ MeV/A'.format(
                beamA=self.beamA,
                beam=self.beam,
                targetA=self.targetA,
                target=self.target,
                beamE=self.beam_energy.value
            )
        else:
            name = f'{self.beam}{self.beamA}{self.target}{self.targetA}E{self.beam_energy.value}'
        return name

    def get_betacms(self):

        if not self._initialize_reaction:
            raise ValueError('Choose the reaction before calculating the betacms.')

        beam_ke = self.beam_energy * self.beamA
        beam_energy_tot = beam_ke + self.beam_mass
        mom_beam = np.sqrt(beam_ke**2 + 2*beam_ke * self.beam_mass)

        gamma = beam_energy_tot/self.beam_mass
        return mom_beam/(gamma*self.beam_mass + self.target_mass)

    def get_rapidity_beam(self):
        if not self._initialize_reaction:
            raise ValueError('Choose the reaction before calculating the beam rapidity.')

        beam_ke = self.beam_energy * self.beamA
        mom_beam = np.sqrt(beam_ke**2 + 2*beam_ke * self.beam_mass)
        return 0.5 * np.log((beam_ke + self.beam_mass + mom_beam) / (beam_ke + self.beam_mass - mom_beam))


        


        