import matplotlib.pyplot as plt
import numpy as np
import math as m
from scipy import integrate

import sys
sys.path.append('../..')
from freddi import Freddi, State, EvolutionResults

DAY = 86400


from numpy import exp, sin

from lmfit import minimize, Parameters
from scipy import interpolate

_Xi_max = 40
_T_iC = 1e8
Mx = 2e33*9.4
GM = 6.673e-8 * Mx
Mdotout = 0
Cirr = 2.9e-4
kerr =  0.4
Mopt = 2.5e33
period = 1.116*DAY
#rout = 1.7e11

q = Mx/Mopt
def returned(path):
    a = np.genfromtxt(path, names = True)
    x1 = (a['tend']/2 + a['tbegin']/2) - (a['tend'][0]/2 + a['tbegin'][0]/2)
    y1 = a['dotM']*1e18
    yerr1 = [a['dotM']-a['b_dotM'],a['B_dotM']- a['dotM']]
    return [x1*DAY, y1, yerr1]

def residual(params, windrad, path):
    """Function to be minimized

    Keyword arguments:
    params -- Parameters involved in th fitting
    x -- The set of values of the moments of observations
    data -- Source luminosity at these moments
    eps_data -- Observation errors
    windrad -- windradarameter for determining the model used
    """
    lis = returned(path)

    rout = params['rout']
    alpha = params['alpha']
    Mdotin = params['Mdotin']
    x0 = params['x0']

    
    default_kwargs = dict(wind=b'no', Mdotout=Mdotout, Mx=Mx, kerr = kerr, alpha = alpha,
            initialcond=b'quasistat', powerorder=1, opacity=b'OPAL', boundcond = b'Tirr', Thot = 1e4, 
            Cirr = Cirr, time=35*DAY, tau=0.35*DAY, Nx=10000, gridscale=b'linear')
    
    def run(**input_kwargs):
        kwargs = default_kwargs.copy()
        kwargs.update(input_kwargs)
        fr = Freddi(**kwargs)
        return fr
    
    
    if windrad == 'wind and tidal':
        frwT = run(wind=b'__Woods_Which_Shields__' , windparams=[_Xi_max, _T_iC], Mdotin = Mdotin, rout = None, Mopt = Mopt, period = period)
        resultT  = frwT.evolve()

        track = interpolate.splrep(resultT.t + x0*DAY, resultT.Mdot_in)
        model = interpolate.splev(lis[0], track)
         
    elif windrad == 'wind no tidal':
        frw = run(rout = rout, F0=Mdotin*np.sqrt(GM*rout), wind=b'__Woods_Which_Shields__' , windparams=[_Xi_max, _T_iC])
        result  = frw.evolve()

        track = interpolate.splrep(result.t + x0*DAY, result.Mdot_in)
        model = interpolate.splev(lis[0], track)
                   
    elif windrad == 'no wind and tidal':
        fr0T = run(Mdotin = Mdotin, rout = None, Mopt = Mopt, period = period)
        r0T  = fr0T.evolve()

        track = interpolate.splrep(r0T.t + x0*DAY, r0T.Mdot_in)
        model = interpolate.splev(lis[0], track)
    elif windrad == 'no wind no tidal':
        fr0 = run(rout = rout, F0=Mdotin*np.sqrt(GM*rout))
        r0  = fr0.evolve()
        
        track = interpolate.splrep(r0.t + x0*DAY, r0.Mdot_in)
        model = interpolate.splev(lis[0], track)
    else:
        print('Something is wrong')
    
    return (lis[1]-model) / lis[2]

params = Parameters()
params.add('rout', min= 1e11, value = 1.7e11, max = 2e11)
params.add('alpha', min = 0.1, value = 0.8, max = 1.0)
params.add('Mdotin', min = 1e18, value = 1e19, max = 1.5e19)
params.add('x0', min= 1.0, value = 3.5, max = 5.0)
params['x0'].vary = False

def fit(wind, path):
    o = minimize(residual, params, args=(wind, path))
    return o

def Cov(Soul):
    """Determination of the covariance matrix

    Keyword arguments:
    Soul --The result of minimization
    """
    S = Soul.covar.copy()
    i, j = np.indices(S.shape)
    Mind = S / np.sqrt(S[i, i] * S[j, j])
    return Mind


def run(**input_kwargs):
    """Element filler

    Keyword arguments:
    **input_kwargs --Input parameter list
    """
    default_kwargs = dict(wind=b'no', Mdotout=Mdotout, Mx=Mx, kerr = kerr,
            initialcond=b'quasistat', powerorder=1, opacity=b'OPAL',  boundcond = b'Tirr', Thot = 1e4, 
            Cirr = Cirr, time=35*DAY, tau=0.35*DAY, Nx=10000, gridscale=b'linear')
    
    kwargs = default_kwargs.copy()
    kwargs.update(input_kwargs)
    fr = Freddi(**kwargs)
    return fr