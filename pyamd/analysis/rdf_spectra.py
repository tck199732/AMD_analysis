import ROOT
import matplotlib.pyplot as plt
from pyamd import PROJECT_DIR
from pyamd.utilities import ame, dataframe, root6
import pandas as pd
import pathlib
import itertools
import collections
import numpy as np

HistoReader = root6.HistogramReader()
df_helper = dataframe.DataFrameHelper()


ROOT.EnableImplicitMT()
ame_table = ame.AME()


class Cuts:
    """ The cuts here should apply to filtered data of any model. 
    """
    @staticmethod
    def cut_uball_multi(br_name='uball_multi', low=1):
        return f'{br_name} >= {low}'

    @staticmethod
    def cut_hira_multi(br_name='hira_multi', low=1):
        return f'{br_name} >= {low}'
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

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
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

        self.column_names = set()
        # check if the cpp functons are set up
        self._setup_cpp_ame_mass = False
        self._cpp_ame_mass = None

        self.rdf = ROOT.RDataFrame(tree, str(path))
        self.column_names.update(self.rdf.GetColumnNames())

    def _define_mass_number(self, rdf, inplace=True):
        rdf = (rdf
               .Define('uball_A', 'uball_N + uball_Z')
               .Define('hira_A', 'hira_N + hira_Z')
               )
        if inplace:
            self.rdf = rdf
        return rdf

    def _define_mass(self, rdf, inplace=True):
        if not self._setup_cpp_ame_mass:
            self._cpp_ame_mass = self._define_mass_cpp()
            self._setup_cpp_ame_mass = True

        self._define_uball_mass(rdf, br_name='uball_mass', inplace=inplace)
        self._define_hira_mass(rdf, br_name='hira_mass', inplace=inplace)
        return rdf

    def _define_uball_mass(self, rdf, br_name='uball_mass', inplace=True):
        rdf = rdf.Define(br_name, f'{self._cpp_ame_mass}(uball_Z, uball_N)')
        self.column_names.add(br_name)
        if inplace:
            self.rdf = rdf
        return rdf

    def _define_hira_mass(self, rdf, br_name='hira_mass', inplace=True):
        rdf = rdf.Define(br_name, f'{self._cpp_ame_mass}(hira_Z, hira_N)')
        self.column_names.add(br_name)
        if inplace:
            self.rdf = rdf
        return rdf

    def _define_mass_cpp(self, maxZ=4, maxN=4, fcn_name='DefineMass'):
        '''
        a. read particle mass from ame.py 
        b. construct cpp script to feed rdf
        '''
        mass_collection = dict()
        for z, n in itertools.product(range(maxZ+1), range(maxN+1)):
            try:
                mass_collection[(n, z)] = ame_table.get_mass(N=n, Z=z).value
            except:
                continue

        cdata = []
        for i, (nz, m), in enumerate(mass_collection.items()):
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

    def _define_momentum(self, rdf, inplace=True):
        rdf = (rdf
               .Define('uball_pmag', 'sqrt(pow(uball_px, 2.) + pow(uball_py, 2.) + pow(uball_pz, 2.))')
               .Define('hira_pmag', 'sqrt(pow(hira_px, 2.) + pow(hira_py, 2.) + pow(hira_pz, 2.))')
               .Define('uball_kinergy', 'sqrt(pow(uball_pmag, 2.) + pow(uball_mass, 2.)) - uball_mass')
               .Define('hira_kinergy', 'sqrt(pow(hira_pmag, 2.) + pow(hira_mass, 2.)) - hira_mass')
               )
        self.column_names.update({'uball_pmag', 'hira_pmag'})
        if inplace:
            self.rdf = rdf
        return rdf

    def _define_angles(self, rdf, inplace=True):
        rdf = (rdf
               .Define('Hira_theta_rdn', 'Hira_theta * TMath::DegToRad()')
               .Define('Hira_phi_rdn', 'Hira_phi * TMath::DegToRad()')
               )
        self.column_names.update({'Hira_theta_rdn', 'Hira_phi_rdn'})
        if inplace:
            self.rdf = rdf
        return rdf

    """
    def _define_kinergy(self, rdf, inplace=True):
        if not 'Hira_mass' in self.column_names:
            rdf = self._define_mass(rdf, inplace)

        rdf = (rdf
               .Define('Hira_kinergy', 'sqrt(pow(Hira_momentum, 2.) + pow(Hira_mass, 2.)) - Hira_mass')
               .Define('Hira_kinergy_perA', 'Hira_kinergy / Hira_A')
               )

        if inplace:
            self.rdf = rdf

        self.column_names.update({'Hira_kinergy', 'Hira_kinergy_perA'})
        return rdf

    def _define_cms(self, rdf, inplace=True):
        gamma = 1. / np.sqrt(1 - self.betacms**2.)

        rdf = (rdf
               .Define('Hira_pzcms', f'{gamma} * (Hira_pz - {self.betacms} * (Hira_kinergy + Hira_mass))')
               .Define('Hira_momentum_cms', 'sqrt(pow(Hira_pzcms,2.) + pow(Hira_pt, 2.))')
               .Define('Hira_kinergy_cms', 'sqrt(pow(Hira_momentum_cms, 2.) + pow(Hira_mass, 2.)) - Hira_mass')
               )
        if inplace:
            self.rdf = rdf
        self.column_names.update(
            {'Hira_pzcms', 'Hira_momentum_cms', 'Hira_kinergy_cms'})
        return rdf

    def _define_rapidity(self, rdf, inplace=True):

        rdf = (rdf
               .Define('Hira_rapidity', '0.5 * log( (Hira_kinergy + Hira_pz + Hira_mass) / (Hira_kinergy - Hira_pz + Hira_mass) )')
               .Define('Hira_rapidity_cms', '0.5 * log( (Hira_kinergy_cms + Hira_pzcms + Hira_mass) / (Hira_kinergy_cms - Hira_pzcms + Hira_mass) )')
               .Define('Hira_rapidity_normalized', f'Hira_rapidity / {self.beam_rapdity}')
               )
        if inplace:
            self.rdf = rdf
        self.column_names.update(
            {'Hira_rapidity', 'Hira_rapidity_cms', 'Hira_rapidity_normalized'})
        return rdf

    def _define_cuts(
            self,
            rdf,
            badmap_version='V1',
            Nc=5,
            multi=1,
            angular_coverage=(30., 75.),
            inplace=True):

        cuts = [
            HiraCuts.cut_4pi_multiplicity(Nc),
            HiraCuts.cut_multiplicity(multi),
            HiraCuts.cut_kinergy(),
            HiraCuts.cut_coverage(angular_coverage),
            HiraCuts.cut_badmap(badmap_version)
        ]
        cuts = ' && '.join(cuts)
        rdf = rdf.Define('Hira_cut', cuts)
        if inplace:
            self.rdf = rdf
        return rdf
    """


# Some simple analysis class are listed below to demonstrate the usage.

class ImpactParameter_Multiplicity(RDataFrame_AMD):
    """ A class for analyzing impact-parameter vs charged-particle multiplicity detected in microball. Since all the data we need would be `b` and `uball_multi`, we directly inherit the class RDataFrame_AMD and construct 2D histogram.
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
        super().__init__(path, tree, beamA, targetA, energy, beam, target)

    def analyze(self, bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], br_names=['uball_multi', 'b'], hname='h2_multi_b'):
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
        pd.DataFrame
        """
        h2 = (self.rdf
              .Filter(f'{br_names[0]} > 0')
              .Histo2D((hname, '', bins[0], *range[0], bins[1], *range[1]), *br_names)
              ).GetValue()
        return HistoReader.hist2d_to_df(h2)


class ImpactParameter_MassNumber(RDataFrame_AMD):
    """ Study the relation between impact-paramter `b` and nuclei mass number `A`. Normally, we would expect more heavy particle in a less central event where it is more likely to have spectator nuclei.
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
        super().__init__(path, tree, beamA, targetA, energy, beam, target)

    def _define_kinematics(self):
        self._define_mass_number(self.rdf, inplace=True)

    def analyze(self, bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], br_names=['uball_A', 'b'], hname='h2_uballA_b'):
        self._define_kinematics()
        h2 = (self.rdf
              .Histo2D((hname, '', bins[0], *range[0], bins[1], *range[1]), *br_names)
              ).GetValue()
        return HistoReader.hist2d_to_df(h2)


class TransverseMomentum_RapidityLab(RDataFrame_AMD):
    """ Study transverse mometum :math: `P_{t}` vs :math: `\hat{y}_{lab}` rapidity / beam-rapidity. Apply cuts to select particle.
    """

    def __init__(self, path, tree='AMD', beamA=40, targetA=58, energy=140, beam='Ca', target='Ni'):
        super().__init__(path, tree, beamA, targetA, energy, beam, target)

    def cut(self):
        pass

    def _define_kinematics(self):
        pass
    
    def analyze(self):
        pass



if __name__ == '__main__':
    rdf = ImpactParameter_MassNumber(
        '/data/kin/amd/feb2022/b10fm/filtered/Ca40Ni58E140_SkM_table3.root')
    df = rdf.analyze()

    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots()
    ax.hist2d(df.x, df.y, weights=df['z'], bins=[25, 100], range=[[-0.5, 24.5], [0., 10.]], cmap='jet', norm=LogNorm())
    fig.savefig('h2_uballA_b.png')
    