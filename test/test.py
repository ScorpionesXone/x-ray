from numpy.testing import assert_allclose, assert_array_equal

import matplotlib.pyplot as plt
import numpy as np
import math as m
import random
import sys
sys.path.append('../..')
from freddi import Freddi, State, EvolutionResults

from numpy import exp

from lmfit import minimize, Parameters
from scipy import interpolate

import doctest
import unittest

from module import *
import module


class MyTestCase(unittest.TestCase):
    routI = 1.25e11
    MdotinI = 1.2e19
    alphaI = 0.75

    default_kwargs = dict(wind=b'no', Mdotout=Mdotout, Mx=Mx, kerr = kerr,
                initialcond=b'quasistat', powerorder=1, opacity=b'OPAL',  boundcond = b'Tirr', Thot = 1e4, 
                Cirr = Cirr, time=35*DAY, tau=0.35*DAY, Nx=10000, gridscale=b'linear')
        
    life = run(wind=b'__Woods_Which_Shields__', windparams=[module._Xi_max, module._T_iC], F0=MdotinI*np.sqrt(GM*routI), alpha = alphaI, rout = routI)
    Sin = life.evolve()


    random.seed(888)

    x1 = Sin.t
    y1 = Sin.Mdot_in

    x = np.array(x1)
    y = np.array(y1) 


    err = np.array([random.random() for i in range(len(y))])
    error = err+4e17

    #yB = y + err*4e17
    #yb = y - err*4e17


    #with open('test.tsv', "w") as file_result:
    #    file_result.write('#tbegin' + ' ' + 'tend' + '\t' + 'dotM' + '\t' + 'b_dotM' + '\t' + 'B_dotM' +  '\n')
    #    for i in range(0, len(y)):
    #        file_result.write(str(x[i]) + '\t' + str(x[i]) + '\t' + str(y[i]) + '\t' + str(yb[i]) + '\t' + str(yB[i]) +  '\n')


    params = Parameters()
    params.add('rout', min= 1e11, value = 1.7e11, max = 4e11)
    params.add('alpha', min = 0.1, value = 0.4, max = 1.6)
    params.add('Mdotin', min = 1e18, value = 1.4e19, max = 1e20)
    params.add('x0', min= 1.0, value = 3.5, max = 5.0)
    params['x0'].vary = False

    Burst = fitnp(params, x, y, error, 'wind no tidal')


    Dalpha = Burst.params['alpha'].value
    DMdotin = Burst.params['Mdotin'].value
    Drout = Burst.params['rout'].value

    print(Dalpha, DMdotin, Drout)

    def test_alpha(self):
        assert_allclose(self.Dalpha, self.alphaI, rtol=0.3)

    def test_rout(self):
        assert_allclose(self.Drout, self.routI, rtol=0.3)
    
    def test_Mdout(self):
        assert_allclose(self.DMdotin, self.MdotinI, rtol=0.3)

if __name__ == '__main__':
    unittest.main()
