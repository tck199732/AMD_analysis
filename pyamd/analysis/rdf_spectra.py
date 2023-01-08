import ROOT
import matplotlib.pyplot as plt
from pyamd import PROJECT_DIR
from pyamd.utilities import ame, dataframe, root6
from pyamd.e15190 import e15190
import pandas as pd
import pathlib
import itertools
import collections
import numpy as np
from typing import Literal

HistoReader = root6.HistogramReader()
df_helper = dataframe.DataFrameHelper()


ROOT.EnableImplicitMT()
ame_table = ame.AME()

def declare_mass_function(maxZ=25, maxN=25, fcn_name='DefineMass'):
    """ A global function for declaring a C++ function which takes `Z` and `N` and return the nuclei mass, which is obtained from the AME table. 
    Parameters
    ----------
    maxZ, maxN : int
        maximum `Z` and `N` to be included in the C++ function. Default value is 25. This value should be large enough to describe data from peripheral collision. 
    fcn_name : str
        The C++ function name which is then used in RDF.Define()
    Return
    ------
    str : fcn_name
    """
    mass_collection = dict()
    for z, n in itertools.product(range(maxZ+1), range(maxN+1)):
        try:
            mass_collection[(n, z)] = ame_table.get_mass(N=n, Z=z).value
        except:
            continue
    cdata = []
    for nz, m in mass_collection.items():
        n, z = nz
        s = f'{{ {{ {n}, {z} }}, {m} }}'
        cdata.append(s)
    cdata = ' ,\n\t'.join(cdata)

    script = f'''
    std::map<std::pair<int, int>, double> AME_MASS_TABLE = {{
        {cdata}
    }};
    ROOT::RVec<double> {fcn_name}(const ROOT::RVec<int>& arrZ, const ROOT::RVec<int>& arrN) {{
        
        ROOT::RVec<double> mass; 
        double m = 1e5;
        for (unsigned int i = 0; i < arrZ.size(); i++) {{
            m = AME_MASS_TABLE[{{arrN[i], arrZ[i]}}];
            mass.push_back(m);
        }}
        return mass;
    }}
    '''
    mangled_name = ROOT.gInterpreter.GetMangledNameWithPrototype(
        0, fcn_name, 'const ROOT::RVec<int>&, const ROOT::RVec<int>&')

    if mangled_name == '':
        ROOT.gInterpreter.Declare(script)

    return fcn_name

class Expr:
    """ A class for constructing expression to be fed into rdf.Define()
    """
    @staticmethod
    def mass_number(Z_name, N_name):
        return f'{Z_name} + {N_name}'

    @staticmethod
    def mass(fcn_name, Z_name, N_name):
        return f'{fcn_name}({Z_name}, {N_name})'
    
    @staticmethod
    def transverse_momentum(px_name, py_name):
        return f'sqrt(pow({px_name}, 2.) + pow({py_name}, 2.))'
    
    @staticmethod
    def momentum(px_name, py_name, pz_name):
        return f'sqrt(pow({px_name}, 2.) + pow({py_name}, 2.) + pow({pz_name}, 2.))'

    @staticmethod
    def kinergy(pmag_name, mass_name):
        return f'sqrt(pow({pmag_name}, 2.) + pow({mass_name}, 2.)) - {mass_name}'
    
    @staticmethod
    def rapidity(kinergy_name, pz_name, mass_name):
        return f'0.5 * log( ({kinergy_name} + {pz_name} + {mass_name}) / ({kinergy_name} - {pz_name} + {mass_name}) )'

    @staticmethod
    def theta_deg(pmag_trans_name, pz_name):
        return f'atan2({pmag_trans_name}, {pz_name}) * TMath::RadToDeg()'
    
    @staticmethod
    def phi_deg(px_name, py_name):
        return f'atan2({py_name}, {px_name}) * TMath::RadToDeg()'
    
    @staticmethod
    def boost_to_lab(beta, pz_name, kinergy_name, mass_name):
        gamma = 1. / np.sqrt(1 - beta**2.)
        return f'{gamma} * ({pz_name} + {beta} * ({kinergy_name} + {mass_name}))'

    @staticmethod
    def boost_to_cms(beta, pz_name, kinergy_name, mass_name):
        gamma = 1. / np.sqrt(1 - beta**2.)
        return f'{gamma} * ({pz_name} - {beta} * ({kinergy_name} + {mass_name}))'

class Cuts:
    """ The cuts here should apply to filtered data of any model. 
    """
    @staticmethod
    def uball_multi(br_name='uball_multi', low=1):
        return f'{br_name} >= {low}'

    @staticmethod
    def hira_multi(br_name='hira_multi', low=1):
        return f'{br_name} >= {low}'

    @staticmethod
    def impact_parameter(br_name='b', low=0., high=10.):
        return f'{br_name} >= {low} && {br_name} < {high}'

    """
    @staticmethod
    def cut_kinergy():
        ekin = {
            'p': (20.0, 198.0),
            'd': (15.0, 263. / 2),
            't': (12.0, 312. / 3),
            '3He': (20.0, 200.0),
            '4He': (18.0, 200.0),
        }
        n = 'hira_N
        z = 'Hira_Z'
        e = 'Hira_kinergy_perA'

        cuts = ' || '.join([
            f'({a} == 1 && {z} == 1 && {e} >= {ekin["p"][0]} && {e} <= {ekin["p"][1]})',
            f'({a} == 2 && {z} == 1 && {e} >= {ekin["d"][0]} && {e} <= {ekin["d"][1]})',
            f'({a} == 3 && {z} == 1 && {e} >= {ekin["t"][0]} && {e} <= {ekin["t"][1]})',
            f'({a} == 1 && {z} == 2 && {e} >= {ekin["3He"][0]} && {e} <= {ekin["3He"][1]})',
            f'({a} == 2 && {z} == 2 && {e} >= {ekin["4He"][0]} && {e} <= {ekin["4He"][1]})',
        ])

        return f'({cuts})'
    """

    @staticmethod
    def select_uball_particle(br_name=['uball_Z', 'uball_N'], values=[1, 0]):
        return '&&'.join([f'{name} == {value}' for name, value in zip(br_name, values)])

    @staticmethod
    def select_hira_particle(br_name=['hira_Z', 'hira_N'], values=[1, 0]):
        return '&&'.join([f'{name} == {value}' for name, value in zip(br_name, values)])


class RDataFrame_AMD:
    """ This class only works on the experimentally-filtered data. In these data sets, there are branches for microball and hira10 detector. The main branches are `N`, `Z`, `px`, `py`, `pz`. Except for the impact parameter, all branches have prefix of the detector names. For instance,
    Branches
    --------
    `b`, `uball_multi`, `hira_px`, ..., etc
    Note
    ----
    Note that the true impact parameter is not the same as the :math: `$\hat{b}$` stored in experiemental data which are determined from `uball_multiplicity`. 
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni', **kwargs):
        """ Initialize the RDataFrame. Default reaction in `Ca40Ni58E140`.
        Parameters
        ----------
        path : str or pathlib.Path
            path of data
        beamA, targetA : int 
            beam / target nuclei mass number
        energy : int 
            beam kinetic energy per nucleon
        beam, target : str
            symbol of beam / target nuclei
        """
        self.beamA, self.targetA, self.energy = beam, target, energy
        self.beam, self.target = beam, target
        self.reaction = f'{beam}{beamA}{target}{targetA}E{energy}'

        self.simulation_detail = dict()
        self.simulation_detail.update(kwargs)

        self.column_names = set()

        self.rdf = ROOT.RDataFrame(tree, str(path))
        self.column_names.update(self.rdf.GetColumnNames())

    def _define_template(self, rdf, rules, inplace=True):
        for br, rule in rules.items():
            if not br in self.column_names:
                rdf = rdf.Define(br, rule)
                self.column_names.add(br)
            else:
                rdf= rdf.Redefine(br, rule)

        if inplace:
            self.rdf = rdf
        return rdf

    def _filter_template(self, rdf, rules, inplace=True):
        rdf = rdf.Filter('&&'.join(rules))
        if inplace:
            self.rdf = rdf
        return rdf

    def define_mass_number(self, rdf, br_name='uball_A', Z_name='uball_Z', N_name='uball_N', inplace=True):
        rules = {
            br_name : Expr.mass_number(Z_name, N_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_mass(self, rdf, fcn_name='DefineMass', br_name='uball_mass', N_name='uball_N', Z_name='uball_Z', inplace=True):
        
        declare_mass_function(fcn_name)
        rules = {
            br_name : Expr.mass(fcn_name, Z_name, N_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)
    
    def define_transverse_momentum(self, rdf, br_name='uball_pmag_trans', px_name='uball_px', py_name='uball_py', inplace=True):
        rules = {
            br_name : Expr.transverse_momentum(px_name, py_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_momentum(self, rdf, br_name='uball_pmag', px_name='uball_px', py_name='uball_py', pz_name='uball_pz', inplace=True):
        rules = {
            br_name : Expr.momentum(px_name, py_name, pz_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_kinergy(self, rdf, br_name='uball_kinergy', pmag_name='uball_pmag', mass_name='mass_name', inplace=True):
        rules = {
            br_name : Expr.kinergy(pmag_name, mass_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_rapidity(self, rdf, br_name='uball_rapidity', kinergy_name='uball_kinergy', mass_name='uball_mass', pz_name='uball_pz', inplace=True):
        rules = {
            br_name : Expr.rapidity(kinergy_name, pz_name, mass_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_theta_deg(self, rdf, br_name='uball_theta_deg', pmag_trans_name='uball_pmag_trans', pz_name='uball_pz', inplace=True):
        rules = {
            br_name : Expr.theta_deg(pmag_trans_name, pz_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_phi_deg(self, rdf, br_name='uball_phi_deg', px_name='uball_px', py_name='uball_py', inplace=True):
        rules = {
            br_name : Expr.phi_deg(px_name, py_name)
        }
        update_rules = {
            br_name : f'{br_name} > 0. ? {br_name} : {br_name} + 360.',
        }
        rules.update(update_rules)
        return self._define_template(rdf, rules, inplace=inplace)

    def define_pz_lab(self, rdf, br_name='uball_pz_lab', pz_name='uball_pz', kinergy_name='kinergy_name', mass_name='uball_mass', beta=None, inplace=True):
        if beta is None:
            beta = e15190.CollisionReaction.get_betacms(self.reaction)
        rules = {
            br_name : Expr.boost_to_lab(beta, pz_name, kinergy_name, mass_name)
        }
        return self._define_template(rdf, rules, inplace=inplace)

    def define_normalized_rapidity_lab(self, rdf, br_name='uball_rapidity_lab_normed', rapidity_name='uball_rapidity_lab', beam_rapidity=None, inplace=True):
        if beam_rapidity is None:
            beam_rapidity = e15190.CollisionReaction.get_rapidity_beam(self.reaction)
        rules = {
            br_name : f'{rapidity_name} / {beam_rapidity}',
        }
        return self._define_template(rdf, rules, inplace=inplace)

# Some simple analysis class are listed below to demonstrate the usage.
class ImpactParameter_Multiplicity(RDataFrame_AMD):
    """ A class for analyzing impact-parameter vs charged-particle multiplicity detected in microball. Since all the data we need would be `b` and `uball_multi`, we directly inherit the class RDataFrame_AMD and construct 2D histogram.
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
        super().__init__(path, tree, beamA, targetA, energy, beam, target)

    def analyze(self, hname='h2_multi_b', bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], br_names=['uball_multi', 'b'], return_type=Literal['ROOT.TH2D', 'pd.DataFrame']='pd.DataFrame'):
        """ Construct 2D histogram Impact-parameter VS multiplicity
        Parameters
        ----------
        hname : str
            name of the TH2D
        bins : sequence of float of size 2
            number of bins in `x`, `y`
        range :  sequence of float of size 2 by 2
            range in `x` and `y`
        br_name : sequence of str of size 2
            branch name of impact-paramter and multiplicity in the experimentally filtered data
        Return
        ------
        pd.DataFrame or ROOT.TH2D
        """
        h2 = (self.rdf
              .Filter(f'{br_names[0]} > 0')
              .Histo2D((hname, '', bins[0], *range[0], bins[1], *range[1]), *br_names)
              ).GetValue()
        if return_type == 'ROOT.TH2D':
            return h2
        elif return_type == 'pd.DataFrame':
            return HistoReader.hist2d_to_df(h2)
        else:
            raise Exception('invalid return_type.')

class TransverseMomentum_RapidityLab(RDataFrame_AMD):
    """ Study transverse mometum :math: `P_{t}` vs :math: `\hat{y}_{lab}` rapidity / beam-rapidity. Apply cuts to select particle.
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
        super().__init__(path, tree, beamA, targetA, energy, beam, target)
        self.results = dict()

    def _select_particles(self, br_name='proton', particle_mass_chage=(1, 1), inplace=False):
        rules = {
            br_name: Cuts.select_hira_particle(values=particle_mass_chage)
        }
        return self._define_template(self.rdf, rules, inplace=inplace)

    def _define_event_cut(self, uball_multi=10, impact_parameter=(0., 10.), inplace=True):
        uball_cut = Cuts.uball_multi(low=uball_multi)
        bcut = Cuts.impact_parameter(
            low=impact_parameter[0], high=impact_parameter[1])
        return self._define_event_cut(self.rdf, [uball_cut, bcut], inplace=inplace)

    def analyze(self, uball_multi=10, impact_parameter=(0., 10.), hname='h2_pta_rapidity_lab_normed', bins=[100, 800], range=[[0., 1.], [0., 800.]]):

        self._define_event_cut(uball_multi, impact_parameter, inplace=True)
        particles = {
            'proton': (1, 1),
            'deuteron': (2, 1),
            'triton': (3, 1),
            'He3': (3, 2),
            'He4': (4, 2),
        }
        beam_rapidity = 1.
        for name, mc in particles.items():
            rdf = self._select_particles(name, *mc, inplace=False)
            h2 = (rdf
                  .Histo((f'{hname}_{name}', '', bins[0], *range[0], bins[1], *range[1]), f'hira_rapidity_lab[{name}] / {beam_rapidity}', f'hira_pmag_trans[{name}]')
                  ).GetValue()
            self.results[name] = h2

    def Output(self, fname=None):
        if fname is None:
            fname = f'{self.reaction}.root'

        f = ROOT.TFile(fname, 'RECREATE')
        f.cd()
        for name, h2 in self.results.items():
            h2.Write()
        f.Write()
        f.Close()


if __name__ == '__main__':
    path = '/data/kin/amd/feb2022/b10fm/filtered/Ca40Ni58E140_SkM_table3.root'
    tree = 'AMD'
    # path = '/data/kin/imqmd-fanurs/filtered/Ca40Ni58E140_cluster.root'
    # tree = 'cluster'
    rdf = ImpactParameter_Multiplicity(
        path=path, tree= tree)
    df = rdf.analyze()

    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots()
    ax.hist2d(df.x, df.y, weights=df['z'], bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], cmap='jet', norm=LogNorm())
    fig.savefig('h2_multi_b.png')

    