# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:35:01 2022
Domaine : 
Sous-Domaine : 
Chapitre :
Appellation : 
@author: Antonin
"""
#-Clean Workspace
# get_ipython().magic('reset -sf')
#
#-Import libraries
import numpy as np
from scipy.integrate import solve_ivp
##
#-Plot figures
import matplotlib.pyplot as plt
#
def TantalOscillator(nf,Np,S,s,H,h0,h,D,tini):
    g = 9.81#m/s**2
    # --- Paramètres ---
    eps = S/s
    print(('epsilon = %0.3f'%eps))
    inveps = eps**-1
    tau = np.sqrt(h/(2*g))
    print(('tau = %0.3f s'%tau))
    D0 = 10**-4*s*np.sqrt(2*g*h)*10**3#L/s
    print(('D0 = %0.3f L/s'%D0))
    D1 = 10**-4*s*np.sqrt(2*g*(h-h0))*10**3#L/s
    print(('D1 = %0.3f L/s'%D1))
    # 
    ratioh = h0/h#h0/h
    etalim = np.sqrt(1-ratioh)
    print('D1/D0 = %0.3f'%etalim)
    # 
    # --- ETUDE DE LA DUREE DE REMPLISSAGE ET DE VIDANGE ---
    # --- Variable = Débit du robinet ---
    #  régime périodique ssi D < D1 <=> eta < etalim
    eta = np.linspace(0.010,0.999999*etalim,Np)#q.epsilon = D/D0
    # 
    # --- Evolution de la durée de la vidange ---
    # w_v = t_v / tau
    # w_v = 2*epsilon*(1-etalim+eta*ln((1-eta)/(etalim-eta)))
    w_v = 2*eps*(1-etalim + eta*np.log((1-eta)/(etalim-eta)))
    # --- Evolution de la durée de remplissage ---
    w_r = (1-etalim**2)*eps*(eta**-1)
    # Plot
    # nf = nf+1
    # plt.figure(nf)
    # plt.semilogy(eta,w_v)
    # plt.semilogy(eta, w_r)
    # plt.semilogy(eta,2*eps*(np.exp(eta) - etalim),'--', color='green')
    # plt.semilogy(eta,(1-etalim**2)*eps*(1-np.exp(-eta))**-1,'--',color='green')
    # plt.xlabel('$D/D_0<D_1/D_0$')
    # plt.xlim([0,etalim])
    # plt.ylabel('$\\tau_v/\\tau, \\tau_r/\\tau$')
    # plt.title(('Evolution de $\\tau_v$ et $\\tau_r$ en fonction du débit $D$ \n $\\epsilon=$ %0.2f'%eps + ' et ' + '$D_1/D_0=$%0.2f'%etalim))
    # plt.legend((('$\\tau_v/\\tau$'),('$\\tau_r/\\tau$'),('évolutions exponentielles')))
    # plt.grid(True,'major','both')
    # 
    """
    Description d'une oscillation de relaxation
    """
    # --- Paramètres ---
    r =  D/D0 
    print('D/D0 =%0.3f'%r)
    q = inveps * r
    # 
    """
    Evolution de la hauteur de la surface libre (phase de remplissage)
    """
    # -
    wini = tini/tau
    # -
    tr = h0*(S*10**-4)/(D*10**-3)
    print("tr = %.2fs"%tr)
    wr = tr/tau
    # -
    ww = np.linspace(wini,wr,Np)
    # -
    coeff = inveps*r
    theta_r = etalim**2 + coeff*ww
    # 
    # - Plot
    # nf = nf+1
    # plt.figure(nf)
    # plt.plot(ww,theta_r)
    # plt.xlabel('$t/\\tau$')
    # plt.ylabel('$z/h$')
    # plt.title('Evolution de la hauteur de la surface libre (phase de remplissage)')
    # plt.grid(True,'major','both')
    # 
    """
    Evolution de la hauteur de la surface libre (phase de vidange)
    """
    # --- Définition de la fonction f pour l'équadiff ---
    # 
    def f(w,theta):
        return q-inveps*np.sqrt(theta)
    # 
    # --- Intervalle de temps
    # 
    tv = 2*eps*tau*(1 - etalim + r * np.log((1 - r)/(etalim - r)))
    print('tv = %0.2fs'%tv)
    # 
    t0 = tr
    tend = tr+tv
    # 
    w0 = t0/tau
    wend = tend/tau
    # 
    w_span = [w0,wend]
    w = np.linspace(w0,wend,Np)
    # 
    # --- Conditions initiales ---
    # 
    Vini = np.array([1])#theta0
    # 
    # --- Résolution numérique ---
    # 
    V = solve_ivp(f, w_span, Vini, method='RK45', t_eval=w, dense_output=False, events=None, vectorized=False, args=None)
    theta_v = V.y[0]
    #
    # --- Résolution approchée ---
    # 
    theta_v_approx = (1+0.5*inveps*(w0-w))**2
    theta_v_approx[theta_v_approx<etalim**2]= etalim**2
    tvapprox = 2*eps*tau*(1-etalim)
    print(('tvapprox = %0.2fs'%tvapprox))
    # 
    # --- Comportement linéaire pour la comparaison ---
    # 
    theta_v_lin = 1 + 0.5*inveps*(w0-w)
    theta_v_lin[theta_v_lin<etalim**2] = etalim**2
    # 
    # --- Plot ---
    # 
    # nf = nf+1
    # plt.figure(nf)
    # plt.plot(w, theta_v, color='purple')
    # plt.xlim([w0, wend])
    # plt.plot(w, theta_v_approx, '--',color='black')
    # plt.plot(w, theta_v_lin, ':',color='blue')
    # plt.hlines(etalim**2,w0,wend,colors='black')
    # plt.xlabel('$t/\\tau$')
    # plt.ylabel('$z/h$')
    # plt.title('Evolution de la hauteur de la surface libre (phase de vidange)')
    # plt.legend((('Solution numérique exacte $D/D_0=$%.2f'%r),'Solution analytique approchée $D/D_0 \ll 1$ (parabole)','Evolution linéaire'))
    # plt.grid(True,'major','both')
    # 
    # --- Etude de l'écart entre la solution approchée ---
    # 
    vecr = eta
    # 
    tvappOvertv = ( 1 + (vecr/(1-etalim)) * np.log((1-vecr)/(etalim-vecr)) )**-1
    # 
    # nf = nf+1
    # plt.figure(nf)
    # plt.title('Ecart sur la durée de vidange entre la solution approchée $D/D_0 \ll 1$ \n et la solution réelle $D/D_0$ qcq')
    # plt.plot(vecr,100*(1-tvappOvertv))
    # plt.plot(vecr,100*vecr,'--',color='black')
    # plt.grid(True,'major','both')
    # plt.xlim([0,etalim])
    # plt.ylim([0,100])
    # plt.xlabel('$D/D_0$')
    # plt.ylabel('$1-\\tau_{v,app}/\\tau_v$[%]')
    # plt.legend((('$1-\\tau_{v,app}/\\tau_v$'),('Evolution linéaire')))
    # 
    """
    Oscillation de relaxation
    """
    # 
    T = tr + tv
    print("Période des Osc. de Rel. : T = %.2fs"%T)
    # 
    theta_T = np.concatenate((theta_r,theta_v))
    w_T = np.concatenate((ww,w))
    # 
    # nf = nf+1
    # plt.figure(nf)
    # plt.title('Une oscillation de relaxation')
    # plt.plot(w_T,theta_T)
    # plt.grid(True,'major','both')
    # plt.xlabel('$t/\\tau$')
    # plt.ylabel('$z/h$')
    # plt.xlim([0,T/tau])
    # plt.legend((['$\\epsilon=$%0.2f'%eps+ ' et ' + '$D/D_0$=%0.2f'%r]))
    # 
    """
    Portrait de phase
    """
    # nf = nf+1
    theta_v_pt = q - inveps*np.sqrt(theta_v)
    theta_r_pt = q * np.ones(np.size(theta_r))
    theta_T_pt = np.concatenate((theta_r_pt,theta_v_pt))
    # 
    # plt.figure(nf)
    # plt.title("Portrait de phase de l'oscillation")
    # plt.plot(theta_T,theta_T_pt)
    # plt.xlabel('$\\theta := z/h$')
    # plt.ylabel('$\\mathrm{d}_w \\theta \quad w = t/\\tau$')
    # plt.grid(True,'major','both')
    # 
    return w_T, theta_T, theta_T_pt, T, tr, etalim, r, eps, tau