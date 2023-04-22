import pathlib
import re
import numpy as np
from astropy import units
from collections import defaultdict
from typing import Literal

from pyamd import PROJECT_DIR
from pyamd.utilities import ame
ame_table = ame.AME()

class Particle:
    ALIAS = {
        'n': 'n1',
        'p': 'H1',
        'd': 'H2',
        't': 'H3',
        '3He': 'He3',
        '4He': 'He4',
        'neutron' : 'n1',
        'proton': 'H1',
        'deuteron': 'H2',
        'triton': 'H3',
        'Helium3': 'He3',
        'Helium4': 'He4',
        'alpha' : 'He4',
        'Alpha' : 'He4',
    }
    def __init__(self, name):
        if name in self.ALIAS:
            name = self.ALIAS[name]
        self.name = name
        self.mass = ame_table.get_mass(name)
        self.N, self.Z = ame_table.get_NZ(symbol=self.name)
    
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


class CollisionReaction:
    @staticmethod
    def dissemble_reaction(reaction, dtype:Literal['dict', 'list']='dict'):
        beamA, targetA, beam_energy = [
            int(m) for m in re.compile('[0-9]+').findall(reaction)]
        beam, target = [
            m for m in re.compile('[A-Za-z]{2}').findall(reaction)]
        
        if dtype == 'dict':
            return {
                'beam' : beam,
                'target' : target,
                'beamA' : beamA,
                'targetA' : targetA,
                'beam_energy' : beam_energy,
            } 
        elif dtype == 'list':
            return beam, target, beamA, targetA, beam_energy
        else:
            raise Exception('`dtype` must be `dict` or `list`.')

    @staticmethod
    def get_betacms(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame_table.get_mass(d['beam'] + str(d['beamA'])) 
        target_mass = ame_table.get_mass(d['target'] + str(d['targetA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        beam_energy_tot = beam_ke + beam_mass
        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)
        gamma = beam_energy_tot / beam_mass
        return (mom_beam / (gamma * beam_mass + target_mass)).value

    @staticmethod
    def get_rapidity_beam(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame_table.get_mass(d['beam'] + str(d['beamA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)
        return 0.5 * np.log((beam_ke + beam_mass + mom_beam) / (beam_ke + beam_mass - mom_beam)).value

    @staticmethod
    def get_rapidity_beam_and_target(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame_table.get_mass(d['beam'] + str(d['beamA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)

        target_mass = ame_table.get_mass(d['target'] + str(d['targetA']))
        
        return (0.5 * np.log((beam_ke + beam_mass + target_mass + mom_beam) / (beam_ke + beam_mass + target_mass - mom_beam))).value

        


        