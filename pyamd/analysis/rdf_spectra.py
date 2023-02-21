import ROOT
import pathlib
import itertools
import numpy as np
import pandas as pd
from typing import Literal

from pyamd.utilities import ame, root6
from pyamd.e15190 import e15190

ROOT.EnableImplicitMT()
HistoReader = root6.HistogramReader()
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
    def magnitude(*args):
        return 'sqrt(' + '+'.join([
            f'pow({name}, 2.)' for name in args
        ]) + ')'

    @staticmethod
    def kinergy(pmag_name, mass_name):
        return f'sqrt(pow({pmag_name}, 2.) + pow({mass_name}, 2.)) - {mass_name}'

    @staticmethod
    def rapidity(kinergy_name, pz_name, mass_name):
        return f'0.5 * log( ({kinergy_name} + {pz_name} + {mass_name}) / ({kinergy_name} - {pz_name} + {mass_name}) )'

    @staticmethod
    def theta_deg(ptrans_name, pz_name):
        return f'atan2({ptrans_name}, {pz_name}) * TMath::RadToDeg()'

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


class RDataFrame_AMD:
    """ This class works on both experimental-filtered data and default data. The only difference is in the momenta data which is stored in the unit of `MeV per nucleon` in the original root file and `MeV` in the filtered file. In filtered-data, branches are named according to the detector microball and hira10 (might add VetoWall and NeutronWall in the future). The main branches are `N`, `Z`, `px`, `py`, `pz`. Except for the impact parameter, all branches have prefix of the detector names. For instance,
    Branches (filtered data)
    ------------------------
    `b`, `uball_multi`, `hira_px`, ..., etc
    Branches (default data)
    -----------------------
    `multi`, `b`, `px`, `py`, `pz`, `N`, `Z`
    Note
    ----
    Note that the true impact parameter is not the same as the :math: `$\hat{b}$` stored in experiemental data which are determined from `uball_multiplicity`. 
    """

    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140'):
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

        self.reaction = reaction

        if isinstance(path, str) or isinstance(path, pathlib.Path):
            self.rdf = ROOT.RDataFrame(tree, str(path))
        elif isinstance(path, list):
            self.rdf = ROOT.RDataFrame(tree, list(map(str, path)))

        self.column_names = set(self.rdf.GetColumnNames())

    def _define_template(self, rules, forced=False, inplace=True):
        for br, rule in rules.items():
            if not br in self.column_names:
                rdf = self.rdf.Define(br, rule)
                self.column_names.add(br)
            elif forced:
                rdf = self.rdf.Redefine(br, rule)
            else:
                return self.rdf
        if inplace:
            self.rdf = rdf
        return rdf

    def _filter_template(self, *args, inplace=True):
        rdf = self.rdf.Filter('&&'.join(args))
        if inplace:
            self.rdf = rdf
        return rdf

    def define_mass_number(self, br_name='A', Z_name='Z', N_name='N', inplace=True):
        rules = {
            br_name: Expr.mass_number(Z_name, N_name)
        }
        return self._define_template(rules, inplace=inplace)

    def redefine_quantity_for_nucleon(self, *args, A_name='A', operator:Literal['*', '/', '+', '-']='*', inplace=True):
        rules = {
            name : f'{name} {operator} {A_name}' for name in args
        }
        return self._define_template(rules=rules, inplace=inplace)

    def define_mass(self, fcn_name='DefineMass', br_name='mass', N_name='N', Z_name='Z', inplace=True):

        declare_mass_function(fcn_name=fcn_name)
        rules = {
            br_name: Expr.mass(fcn_name, Z_name, N_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_transverse_momentum(self, br_name='ptrans', px_name='px', py_name='py', inplace=True):
        rules = {
            br_name: Expr.magnitude(px_name, py_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_momentum(self, br_name='pmag', px_name='px', py_name='py', pz_name='pz', inplace=True):
        rules = {
            br_name: Expr.magnitude(px_name, py_name, pz_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_kinergy(self, br_name='kinergy', pmag_name='pmag', mass_name='mass', inplace=True):
        rules = {
            br_name: Expr.kinergy(pmag_name, mass_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_rapidity(self, br_name='rapidity', kinergy_name='kinergy', mass_name='mass', pz_name='pz', inplace=True):
        rules = {
            br_name: Expr.rapidity(kinergy_name, pz_name, mass_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_theta_deg(self, br_name='theta_deg', ptrans_name='ptrans', pz_name='pz', inplace=True):
        rules = {
            br_name: Expr.theta_deg(ptrans_name, pz_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_phi_deg(self, br_name='phi_deg', px_name='px', py_name='py', forced=True, inplace=True):
        rules = {
            br_name: Expr.phi_deg(px_name, py_name)
        }
        update_rules = {
            br_name: f'{br_name} > 0. ? {br_name} : {br_name} + 360.',
        }
        rules.update(update_rules)
        return self._define_template(rules, forced=forced, inplace=inplace)

    def define_pz_lab(self, br_name='pz_lab', pz_name='pz', kinergy_name='kinergy', mass_name='mass', beta=None, inplace=True):
        if beta is None:
            beta = e15190.CollisionReaction.get_betacms(self.reaction)
        rules = {
            br_name: Expr.boost_to_lab(beta, pz_name, kinergy_name, mass_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_normalized_rapidity(self, br_name='rapidity_lab_normed', rapidity_name='rapidity_lab', beam_rapidity=None, inplace=True):
        if beam_rapidity is None:
            beam_rapidity = e15190.CollisionReaction.get_rapidity_beam(
                self.reaction)
        rules = {
            br_name: f'{rapidity_name} / {beam_rapidity}',
        }
        return self._define_template(rules, inplace=inplace)

    def define_particle(self, br_name, Z_name='Z', N_name='N', inplace=True):
        par = e15190.Particle(br_name)
        rules = {
            br_name: f'{Z_name} == {par.Z} && {N_name} == {par.N}'
        }
        return self._define_template(rules, inplace=inplace)

    def define_position(self, br_name='rmag', x_name='x', y_name='y', z_name='z', inplace=True):
        rules = {
            br_name: Expr.magnitude(x_name, y_name, z_name)
        }
        return self._define_template(rules, inplace=inplace)

    def define_charged_particles(self, br_name='charged_particle', Z_name='Z', inplace=True):
        rules = {
            br_name : f'{Z_name} != 0'
        }
        return self._define_template(rules, inplace=inplace)

"""
From here on are some analysis classes utilizing the above RDataFrame. User should get the instance of the analysis in a seperate script and call the `analyze` method to get a `ROOT.TH2D` object or pandas.dataframe. 
Analysis for filtered data
--------------------------
1. ImpactParameter_Multiplicity : 2D histogram `impact-parameter vs uball charged-multiplicity`
2. TransverseMomentum_RapidityLab : 2D histogram `Pt vs rapidity / beam-rapidity (lab)`

Analysis for pure data
----------------------
1. Stopping : 1D histogram of total Et / Ez of each event. 
2. TransverseMomentum_RapidityLab : 2D histogram `Pt vs rapidity / beam-rapidity (lab)`
3. SpaceTime : 2D histogram `Emission Time vs Momentum`
"""

class ImpactParameter_Multiplicity(RDataFrame_AMD):
    """ A class for analyzing impact-parameter vs charged-particle multiplicity detected in microball. Since all the data we need would be `b` and `uball_multi`, we directly inherit the class RDataFrame_AMD and construct 2D histogram.
    """

    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140'):
        super().__init__(path, tree, reaction)

    def analyze(self, hname='h2_multi_b', bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], br_names=['uball_multi', 'b'], return_type: Literal['ROOT.TH2D', 'pd.DataFrame'] = 'ROOT.TH2D'):
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

        # 1 primary event -> 10 secondary events
        weight = 0.1

        # nevents = self.rdf.Filter(f'{br_names[0]} > 0').Count().GetValue() * weight
        # bug fix, should normalize by pure simulation event
        nevents = self.rdf.Count().GetValue() * weight
        h2 = (self.rdf
              .Filter(f'{br_names[0]} > 0')
              .Define('weight', f'{weight}')
              .Histo2D((hname, '', bins[0], *range[0], bins[1], *range[1]), *br_names, 'weight')
              ).GetValue()

        h2.Scale(1./nevents)
        if return_type == 'ROOT.TH2D':
            return h2
        elif return_type == 'pd.DataFrame':
            return HistoReader.hist2d_to_df(h2)
        else:
            raise Exception('invalid return_type.')


# this class is only for testing purpose. It is more convenient to analyze in C++ as we need to replace errorbars later.
class TransverseMomentum_RapidityLab(RDataFrame_AMD):
    """ Study transverse mometum :math: `P_{t}` vs :math: `\hat{y}_{lab}` rapidity / beam-rapidity. Same calculations are done for the light clusters `p`, `d`, `t`, `3He` and `4He`. Note that in this analysis the filtered data table3 is read. This means 10 decay events are simulated out of 1 primary event and thus the error calculation is underestimated. For correct error calculation, one can approximate by the error bar obtained from only analyzing 1 decay event for each primary event (See **`${Project_Dir}/analysis/spectra.cpp`**) 
    """

    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140', dtype:Literal['filtered', 'raw']='filtered'):
        super().__init__(path, tree, reaction)
        self.dtype = dtype

    def _define_kinematics(self):
        if self.dtype == 'filtered':
            self._define_kinematics_filtered_data()
        elif self.dtype == 'raw':
            self._define_kinematics_raw_data()
        else:
            raise ValueError('Invalid `dtype`, must be `filtered` or `raw`.')

    def _define_kinematics_filtered_data(self):
        # quantities which are the same lab and cms frame
        self.define_mass_number(
            br_name='hira_A', Z_name='hira_Z', N_name='hira_N')
        self.define_mass(br_name='hira_mass', N_name='hira_N', Z_name='hira_Z')
        self.define_transverse_momentum(
            br_name='hira_ptrans', px_name='hira_px', py_name='hira_py')

        # cms frame quantities
        self.define_momentum(
            br_name='hira_pmag', px_name='hira_px', py_name='hira_py', pz_name='hira_pz')
        self.define_kinergy(br_name='hira_kinergy',
                            pmag_name='hira_pmag', mass_name='hira_mass')

        # lab frame quantities
        self.define_pz_lab(br_name='hira_pz_lab', pz_name='hira_pz',
                           kinergy_name='hira_kinergy', mass_name='hira_mass')
        self.define_momentum(br_name='hira_pmag_lab', px_name='hira_px',
                             py_name='hira_py', pz_name='hira_pz_lab')
        self.define_kinergy(br_name='hira_kinergy_lab',
                            pmag_name='hira_pmag_lab', mass_name='hira_mass')
        self.define_rapidity(br_name='hira_rapidity_lab', kinergy_name='hira_kinergy_lab',
                             mass_name='hira_mass', pz_name='hira_pz_lab')
        self.define_normalized_rapidity(
            br_name='hira_rapidity_lab_normed', rapidity_name='hira_rapidity_lab')


    def _define_kinematics_raw_data(self):
        # quantities which are the same lab and cms frame
        self.define_mass_number(br_name='A', Z_name='Z', N_name='N')
        self.define_mass(br_name='mass', N_name='N', Z_name='Z')

        # change unit to total momenta
        self.redefine_quantity_for_nucleon('px', 'py', 'pz', A_name='A', operator='*')

        self.define_transverse_momentum(br_name='ptrans', px_name='px', py_name='py')

        # cms frame quantities
        self.define_momentum(br_name='pmag', px_name='px', py_name='py', pz_name='pz')
        self.define_kinergy(br_name='kinergy', pmag_name='pmag', mass_name='mass')

        # lab frame quantities
        self.define_pz_lab(br_name='pz_lab', pz_name='pz', kinergy_name='kinergy', mass_name='mass')
        self.define_momentum(br_name='pmag_lab', px_name='px', py_name='py', pz_name='pz_lab')
        self.define_kinergy(br_name='kinergy_lab', pmag_name='pmag_lab', mass_name='mass')
        self.define_rapidity(br_name='rapidity_lab', kinergy_name='kinergy_lab', mass_name='mass', pz_name='pz_lab')
        self.define_normalized_rapidity(br_name='rapidity_lab_normed', rapidity_name='rapidity_lab')
    
    def _filter_event(self):
        if self.dtype == 'filtered':
            self._filter_event_filtered_data()
        elif self.dtype == 'raw':
            self._filter_event_raw_data()

    def _filter_event_raw_data(self, multi=(0, 100), bcut=(0., 10.), inplace=True):
        uball_cut = f'multi >= {multi[0]} && multi <= {multi[1]}'
        bcut = f'b >= {bcut[0]} && b < {bcut[1]}'
        return self._filter_template(uball_cut, bcut, inplace=inplace)


    def _filter_event_filtered_data(self, multi=(10, 25), bcut=(0., 10.), inplace=True):
        uball_cut = f'uball_multi >= {multi[0]} && uball_multi <= {multi[1]}'
        bcut = f'b >= {bcut[0]} && b < {bcut[1]}'
        return self._filter_template(uball_cut, bcut, inplace=inplace)

    def analyze(self, multi=(10, 25), bcut=(0., 10.), hname='h2_pta_rapidity_lab_normed', bins=[100, 800], range=[[0., 1.], [0., 800.]]):
        """ Fill 2D Histograms of Pt/A vs rapidity / beam-rapidity in the lab frame. 
        Paramters
        ---------
        multi : sequence of float of size 2
            range of microball multiplicity
        bcut : sequence of float of size 2
            range of impact paramter. for simulation it is almost always kept unconstrained.
        Return
        ------
        dict of ROOT.TH2D
        """

        self._filter_event(multi, bcut, inplace=True)

        self._define_kinematics()

        particles = [
            'proton',
            'deuteron',
            'triton',
            'Helium3',
            'Helium4',
        ]

        results = dict()
        nevents = self.rdf.Count().GetValue()

        for par in particles:
            self.define_particle(br_name=par, Z_name='hira_Z',
                                 N_name='hira_N', inplace=True)
            A = e15190.Particle(par).Z + e15190.Particle(par).N
            h2 = (self.rdf
                  .Define(f'hira_ptrans_{par}', f'hira_ptrans[{par}]')
                  .Define(f'hira_ptrans_per_A_{par}', f'hira_ptrans_{par} / {A}')
                  .Define(f'hira_rapidity_lab_normed_{par}', f'hira_rapidity_lab_normed[{par}]')
                  .Define(f'weight', '0.1')
                  .Histo2D((f'{hname}_{par}', '', bins[0], *range[0], bins[1], *range[1]), f'hira_rapidity_lab_normed_{par}', f'hira_ptrans_per_A_{par}', 'weight')
                  ).GetValue()

            h2.Scale(10./nevents)
            results[par] = h2

        return results

class Spacetime_Momentum(RDataFrame_AMD):
    """ Study `emission time` and `emission size` vs `Momentum`. The spacetime information can only be extracted from primary data with suffix `table21t.root`. 
    """
    PARTICLES = [
        'proton',
        'deuteron',
        'triton',
        'Helium3',
        'Helium4',
    ]

    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140'):
        super().__init__(path, tree, reaction)

    def _define_kinematics(self):
        """ Be careful here, make sure the unit of `p` is correct (`MeV` or `MeV per nucleon`).
        """
        # quantities which are the same lab and cms frame
        self.define_transverse_momentum(
            br_name='ptrans', px_name='px', py_name='py')

        # cms frame quantities
        self.define_momentum(br_name='pmag', px_name='px',
                             py_name='py', pz_name='pz')
        self.define_position(br_name='rmag', x_name='x',
                             y_name='y', z_name='z')

    def analyze_emission_time(self, hname='h2_time_momentum', bins=[600, 500], range=[[0., 600.], [0., 500.]]):

        self._define_kinematics()
        results = dict()
        nevents = self.rdf.Count().GetValue()

        for par in self.PARTICLES:
            self.define_particle(br_name=par, Z_name='Z',
                                 N_name='N', inplace=True)
            h2 = (self.rdf
                  .Define(f'pmag_{par}', f'pmag[{par}]')
                  .Define(f't_{par}', f't[{par}]')
                  .Define('tcut', f't_{par} > 0.')
                  # note : order of [tcut] is important here
                  .Redefine(f'pmag_{par}', f'pmag_{par}[tcut]')
                  .Redefine(f't_{par}', f't_{par}[tcut]')
                  .Histo2D((f'{hname}_{par}', '', bins[0], *range[0], bins[1], *range[1]), f'pmag_{par}', f't_{par}')
                  ).GetValue()

            h2.Scale(1./nevents)
            results[par] = h2

        return results

    def analyze_spacetime(self, hname='h2_position_momentum', bins=[400, 500], range=[[0., 40.], [0., 500.]]):

        self._define_kinematics()
        results = dict()
        nevents = self.rdf.Count().GetValue()

        for par in self.PARTICLES:
            self.define_particle(br_name=par, Z_name='Z',
                                 N_name='N', inplace=True)

            h2 = (self.rdf
                  .Define(f'rmag_{par}', f'rmag[{par}]')
                  .Define(f't_{par}', f't[{par}]')

                  .Define('tcut', f't_{par} > 0.')
                  # note : order of [tcut] is important here
                  .Redefine(f'rmag_{par}', f'rmag_{par}[tcut]')
                  .Redefine(f't_{par}', f't_{par}[tcut]')

                  .Histo2D((f'{hname}_{par}', '', bins[0], *range[0], bins[1], *range[1]), f'rmag_{par}', f't_{par}')
                  ).GetValue()

            h2.Scale(1./nevents)
            results[par] = h2

        return results

class Stopping(RDataFrame_AMD):
    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140'):
        super().__init__(path, tree, reaction)

    def _define_kinematics(self):
        self.define_mass_number(
            br_name='A', Z_name='Z', N_name='N')

        self.define_mass(br_name='mass', N_name='N', Z_name='Z')
        self.redefine_quantity_for_nucleon('px', 'py', 'pz', A_name='A', operator='*')
        
        self.define_transverse_momentum(br_name='ptrans')
        self.define_kinergy(br_name='etrans', pmag_name='ptrans', mass_name='mass')
        self.define_kinergy(br_name='elong', pmag_name='pz', mass_name='mass')
        self.define_charged_particles(br_name='charged_particle', Z_name='Z')

    def analyze(self, hname='h1_stopping', bins=1000, range=[0.0,2.0], weight=1.):
        
        self._define_kinematics()
        nevents = self.rdf.Count().GetValue()
        h1 = (self.rdf.
            Define('etrans_charged_particle', f'etrans[charged_particle]').
            Define('elong_charged_particle', f'elong[charged_particle]').
            Define('total_etrans_charged_particle', f'ROOT::VecOps::Sum(etrans_charged_particle)').
            Define('total_elong_carged_particle', f'ROOT::VecOps::Sum(elong_charged_particle)').
            Define('stopping', 'total_etrans_charged_particle / total_elong_carged_particle').
            Define('weight', f'{weight}').
            Histo1D((hname, '', bins, *range), 'stopping', 'weight')
        ).GetValue()

        h1.Scale(1. / (nevents * weight))

        return h1

class Rapidity_Distribution(RDataFrame_AMD):
    def __init__(self, path, tree='AMD', reaction='Ca40Ni58E140'):
        super().__init__(path, tree, reaction)
    
    def _define_kinematics(self):
        
        # change unit to total momenta
        self.redefine_quantity_for_nucleon('px', 'py', 'pz', A_name='A', operator='*')

        # quantities which are the same lab and cms frame
        self.define_mass(br_name='mass', N_name='N', Z_name='Z')
        
        # cms frame quantities
        self.define_momentum(
            br_name='pmag', px_name='px', py_name='py', pz_name='pz')
        self.define_kinergy(br_name='kinergy',
                            pmag_name='pmag', mass_name='mass')

        # lab frame quantities
        self.define_pz_lab(br_name='pz_lab', pz_name='pz',
                           kinergy_name='kinergy', mass_name='mass')
        self.define_momentum(br_name='pmag_lab', px_name='px',
                             py_name='py', pz_name='pz_lab')
        self.define_kinergy(br_name='kinergy_lab',
                            pmag_name='pmag_lab', mass_name='mass')
        self.define_rapidity(br_name='rapidity_lab', kinergy_name='kinergy_lab',
                             mass_name='mass', pz_name='pz_lab')

        self.define_normalized_rapidity(
            br_name='rapidity_lab_normed', rapidity_name='rapidity_lab')
    
    def analyze(self, hname='h1_rapidity_lab_normed', bins=400, range=[-2.,2.], weight=1.):
        
        self._define_kinematics()
        nevents = self.rdf.Count().GetValue()

        results = dict()
        particles = [
            'proton',
            'deuteron',
            'triton',
            'Helium3',
            'Helium4',
        ]

        for par in particles:
            self.define_particle(br_name=par, Z_name='Z', N_name='N', inplace=True)
            h1 = (self.rdf.
                Define('weight', f'{weight}').
                Define(f'rapidity_lab_normed_{par}', f'rapidity_lab_normed[{par}]').
                Histo1D((f'{hname}_{par}', '', bins, *range), f'rapidity_lab_normed_{par}', 'weight')
            ).GetValue()

            h1.Scale(1. / (nevents * weight))
            results[par] = h1

        return results