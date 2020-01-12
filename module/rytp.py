import matplotlib.pyplot as plt
import numpy as np
import math as m
from scipy import integrate

import sys
sys.path.append('../..')
pat = sys.argv[1]
from freddi import Freddi, State, EvolutionResults

DAY = 86400


from numpy import exp, sin

from lmfit import minimize, Parameters
from scipy import interpolate
from module import *
import module

def Minimum():
    """Function to be minimize the data"""

    params = Parameters()
    params.add('rout', min= 1e11, value = 1.7e11, max = 2e11)
    params.add('alpha', min = 0.1, value = 0.8, max = 1.0)
    params.add('Mdotin', min = 1e18, value = 1e19, max = 1.5e19)
    params.add('x0', min= 1.0, value = 3.5, max = 5.0)
    params['x0'].vary = False

    out = fit(params, 'wind no tidal', path=pat)
    out.params

    params['rout'].vary = False

    ouT = fit(params, 'wind and tidal', path=pat)
    ouT.params

    ouT0 = fit(params, 'no wind and tidal', path=pat)
    ouT0.params

    lil = []

    lil.append(out.params['alpha'].value)
    
    lil.append(out.params['Mdotin'].value)
    
    lil.append(out.params['rout'].value)

    
    lil.append(ouT.params['alpha'].value)
    
    lil.append(ouT.params['Mdotin'].value)

    
    lil.append(ouT0.params['alpha'].value)
    
    lil.append(ouT0.params['Mdotin'].value)

    return lil

    #alpha0 = out0.params['alpha'].value
    #Mdotin0 = out0.params['Mdotin'].value
    #x00= out0.params['x0'].value

        

def Plotting(rar):
    """Function to be plot the result"""

    #fr0 = run(F0=Mdotin0*np.sqrt(GM*rout), alpha = alpha0, rout = rout)
    #r0  = fr0.evolve()


    #params['alpha'].min =  1.2
    #params['alpha'].value =  1.3
    #params['Mdotin'].min =  1e19
    #params['alpha'].value =  1e19
    #out0 = minimize(residual, params, args=(x1*DAY, y1, yerr1 ,'no wind no tidal'))
    #out0.params

    frT0 = run(alpha = rar[5], Mdotin = rar[6], rout = None, Mopt = Mopt, period = period)
    rT0  = frT0.evolve()

    frw = run(wind=b'__Woods_Which_Shields__', windparams=[module._Xi_max, module._T_iC], F0=rar[1]*np.sqrt(GM*rar[2]), alpha = rar[0], rout = rar[2])
    result = frw.evolve()

    frwT = run(wind=b'__Woods_Which_Shields__', windparams=[module._Xi_max,module. _T_iC], Mdotin = rar[4], alpha = rar[3], rout = None, Mopt = Mopt, period = period)
    resultT = frwT.evolve()


    alphaT1 = round(rar[3], 3)
    alpha0T1 = round(rar[5], 3)

    lis = returned(pat)


    plt.figure(figsize = (10,6))
    plt.title(r'Wind off: $\alpha = $' + str(alphaT1)+r', ' + r'Wind on: $\alpha = $' + str(alpha0T1) + r'; ' + r'Woods Approx Case', fontsize=16)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel(r'$t$, days after peak', fontsize=33)
    plt.ylabel(r'$\dotM$, g/s', fontsize=33)
    plt.errorbar(lis[0]/DAY, lis[1], lis[2], fmt='x', color = 'k', label='Observe')
    plt.plot(rT0.t/DAY + 3.5, rT0.Mdot_in, label='Wind is off')
    plt.plot(resultT.t / DAY + 3.5, resultT.Mdot_in, label='Wind is on')  
    plt.axhline(np.exp(-1)*resultT.Mdot_in[0], ls='-.', color='k', lw=0.5, label='$\mathrm{e}^{-1}$')
    plt.legend()
    plt.grid()
    #plt.savefig('MdotvsTime.pdf', bbox_inches = 'tight')


    plt.figure(figsize = (10,6))
    plt.title(r'Wind off: $\alpha = $' + str(alphaT1)+r', ' + r'Wind on: $\alpha = $' + str(alpha0T1) + r'; ' + r'Woods Approx Case', fontsize=16)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel(r'$t$, days after peak', fontsize=33)
    plt.ylabel(r'$R_{\rm hot}$, cm', fontsize=33) 
    plt.plot(rT0.t[1:]/DAY + 3.5, rT0.last_R[1:], label='Wind is off', color = 'r') 
    plt.errorbar(resultT.t[1:] / DAY + 3.5, resultT.last_R[1:], fmt='*', color = 'k', label='Wind is on') 
    plt.legend()
    plt.grid()
    #plt.savefig('RoutvsTime.pdf', bbox_inches = 'tight')


    plt.figure(figsize = (10,6))
    plt.xlabel(r'$t$, days after peak', fontsize=33)
    plt.ylabel(r'$\dot{M}_{\rm wind}/\dot{M}_{\rm acc}$', fontsize=33)
    plt.plot(resultT.t[1:] / DAY+ 3.5, resultT.Mdot_wind[1:]/resultT.Mdot_in[1:], color = 'g')
    plt.grid()
    #plt.savefig('RelatvsTime.pdf', bbox_inches = 'tight')


    m_P = 1.6726e-24
    k_B = 1.3807e-16
    mu = 0.61
    Ric = ((GM*mu*m_P)/(k_B*module._T_iC))
    SMTH = -resultT.windC[1,:]*GM*GM/(4*m.pi*(resultT.h[1,:])*(resultT.h[1,:])*(resultT.h[1,:]))

    plt.figure(figsize = (10,6))
    plt.xlim(0, resultT.last_R[2]/Ric)
    plt.xlabel(r'$R/R_{\rm IC}$', fontsize=33)
    plt.ylabel(r'$W$, g/(s*cm$^{2}$)', fontsize=33)
    plt.plot(resultT.R[1,:]/Ric, SMTH)
    plt.grid()
    #plt.savefig('WindvsRad.pdf', bbox_inches = 'tight')


def main():

    lol = Minimum()
    Plotting(lol)

if __name__ == '__main__':
    main()
