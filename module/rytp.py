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


def main():

    out = fit('wind no tidal', path=pat)
    out.params

    params['rout'].vary = False

    ouT = fit('wind and tidal', path=pat)
    ouT.params

    ouT0 = fit('no wind and tidal', path=pat)
    ouT0.params


    #params['alpha'].min =  1.2
    #params['alpha'].value =  1.3
    #params['Mdotin'].min =  1e19
    #params['alpha'].value =  1e19
    #out0 = minimize(residual, params, args=(x1*DAY, y1, yerr1 ,'no wind no tidal'))
    #out0.params


        
    print(Cov(out)) 


    alpha = out.params['alpha'].value
    Mdotin = out.params['Mdotin'].value
    rout = out.params['rout'].value

    alphaT = ouT.params['alpha'].value
    MdotinT = ouT.params['Mdotin'].value


    alpha0T = ouT0.params['alpha'].value
    Mdotin0T = ouT0.params['Mdotin'].value

    #alpha0 = out0.params['alpha'].value
    #Mdotin0 = out0.params['Mdotin'].value
    #x00= out0.params['x0'].value

        
    #fr0 = run(F0=Mdotin0*np.sqrt(GM*rout), alpha = alpha0, rout = rout)
    #r0  = fr0.evolve()

    frT0 = run(alpha = alpha0T, Mdotin = Mdotin0T, rout = None, Mopt = Mopt, period = period)
    rT0  = frT0.evolve()

    frw = run(wind=b'__Woods_Which_Shields__', windparams=[module._Xi_max, module._T_iC], F0=Mdotin*np.sqrt(GM*rout), alpha = alpha, rout = rout)
    result = frw.evolve()

    frwT = run(wind=b'__Woods_Which_Shields__', windparams=[module._Xi_max,module. _T_iC], Mdotin = MdotinT, alpha = alphaT, rout = None, Mopt = Mopt, period = period)
    resultT = frwT.evolve()


    alphaT1 = round(alphaT, 3)
    alpha0T1 = round(alpha0T, 3)

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

if __name__ == '__main__':
    main()