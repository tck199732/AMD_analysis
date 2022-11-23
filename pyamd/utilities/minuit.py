

import pathlib
import os
import subprocess
import sys
import numpy as np
from array import array
import ctypes

try:
    import ROOT
except:
    path_root = subprocess.check_output('which root', shell=True)
    path_root = str(path_root, encoding='utf-8').strip()
    path_root = pathlib.Path(str(path_root)).parent.parent
    path_root = pathlib.Path(path_root, 'lib')
    sys.path.append(str(path_root))
    import ROOT


class TMinuit:
    def __init__(self, n, fcn, maxcall=500, tol=0.1):
        self.minuit = ROOT.TMinuit(n)
        self.minuit.SetFCN(fcn)

        self.maxcall = maxcall
        self.tol = tol

        self.arglist = array('d', 10*[0.])
        self.ierflg = ctypes.c_int(0)
        self.arglist[0] = 1
        self.minuit.mnexcm("SET ERR", self.arglist, 1, self.ierflg)

        self.parameters = dict()
        self.errors = dict()
        self.parameter_names = dict()
        self.performance = dict()

    def set_verbose(self, verbose=None):
        if verbose is None:
            verbose = 0
        # -1 : quiet
        self.minuit.SetPrintLevel(verbose)

    def set_parameter(self, n, name, v0, step, range=(0., 0.), fixed=False):
        self.minuit.mnparm(n, name, v0, step, *range, self.ierflg)
        self.parameters[name] = 0.
        self.errors[name] = 0.
        self.parameter_names[n] = name
        if fixed:
            self.minuit.FixParameter(n)

    def fit(self, fitter='MIGRAD', verbose=None):
        self.arglist[0] = self.maxcall
        self.arglist[1] = self.tol
        self.set_verbose(verbose)
        self.minuit.mnexcm(fitter, self.arglist, 2, self.ierflg)

        # save the result here
        # FMIN: the best function value found so far
        # FEDM: the estimated vertical distance remaining to minimum
        # ERRDEF: the value of UP defining parameter uncertainties
        # NPARI: the number of currently variable parameters
        # NPARX: the highest (external) parameter number defined by user
        # ISTAT: a status integer indicating how good is the covariance matrix:
        # 0= not calculated at all
        # 1= approximation only, not accurate
        # 2= full matrix, but forced positive-definite
        # 3= full accurate covariance matrix

        Double = ctypes.c_double
        FMIN, FEDM, ERRDEF = Double(), Double(), Double()
        NPARI, NPARX, ISTAT = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        self.minuit.mnstat(FMIN, FEDM, ERRDEF, NPARI, NPARX, ISTAT)
        NPARI, NPARX, ISTAT = map(lambda x: x.value, [NPARI, NPARX, ISTAT])

        self.performance = {
            'FMIN': FMIN.value,
            'FEDM': FEDM.value,
            'ERRDEF': ERRDEF.value,
            'NPARI': NPARI,
            'NPARX': NPARX,
            'ISTAT': ISTAT,
        }

        # print(self.minuit.GetNumPars() == NPARI)
        for i in range(self.minuit.GetNumPars()):
            par, err = Double(), Double()
            self.minuit.GetParameter(i, par, err)
            name = self.parameter_names[i]
            self.parameters[name] = par.value
            self.errors[name] = err.value

    def get_chisq(self):
        return self.performance['FMIN']

    def get_parameter(self, name):
        return self.parameters[name]

    def get_error(self, name):
        return self.errors[name]

    def get_parameter_name(self, n=0):
        return self.parameter_names[n]

# class IMinuit:
#     def __init__(self, n, fcn):
#         self.minuit = iminuit.Minuit()
#         self.fcn = fcn


# testing
'''
x = array('d', [1.,2.,3.,4.,5.])
y = array('d', [10.,20.,30.,41.,43.])
yerr = array('d', [1.,1.,1.,4.,8.])

def fcn( npar, gin, f, par, iflag):
    chisq = 0.
    for i in range(len(x)):
        delta = (model(x[i], par) - y[i]) / yerr[i]
        chisq += delta**2
    f.value = chisq

def model( x, par ):
    return par[0]*x + par[1]

if __name__ == '__main__':
    gMinuit = TMinuit(2, fcn)
    gMinuit.set_parameter(0, 'a1', 7, 0.1)
    gMinuit.set_parameter(1, 'a2', 0.0, 0.1)
    gMinuit.fit()

'''


m = ROOT.TMinuit(2)


def model(x, par):
    return par[0]*x + par[1]


def fcn(npar, gin, f, par, iflag):
    chisq = 0.
    for i in range(len(x)):
        delta = (model(x[i], par) - y[i]) / yerr[i]
        chisq += delta**2
    f.value = chisq
