# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:40:34 2022
Domaine : 
Sous-Domaine : 
Chapitre :
Appellation : 
@author: Antonin
"""
#-Clean Workspace
get_ipython().magic('reset -sf')
#
#-Import libraries
import numpy as np
from scipy.integrate import solve_ivp
##
#-Plot figures
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.close('all')
# mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.family'] = 'STIXGeneral'
# mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['figure.dpi'] = 110
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.grid.axis'] = 'both'
mpl.rcParams['axes.grid.which'] = 'minor'
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['figure.subplot.hspace']=  0.6
mpl.rcParams['figure.subplot.wspace']=  0.4
mpl.rcParams['lines.linewidth']  = 1.4
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['figure.figsize'] = 1.5*7.4, 1.5*5.8
# mpl.rcParams['text.usetex'] = True
from mpl_toolkits.mplot3d import Axes3D
from LIBvasetantale import TantalOscillator
#
"""
Début
"""
# 
nf = 0#pour numeros figures
#
Np = 10000#Nb de pts pour les plots
# 
#
# --- Constantes physiques du problème ---
# z0 = h 
# tau = sqrt(h/2g)
# q = D*tau/(S*h) = (D/D0) * (1/eps)
# eps = S/s >>1
# r = 1/eps <<1
# -----
# g = 9.81#m/s**2
# 
# #############################################################################

S = 80#cm**2
s = 1#cm**2
# -----
H = 1/2#m
h0 = 3/5*H#m
h = 1.2#m
# 
D = 0.12#LL/s
# 
tini = 0
# 
# #############################################################################
# 
# --- Oscillation de relaxation
# 
w_T, theta_T, theta_T_pt, T, tr, etalim, r, eps, tau = TantalOscillator(nf,Np,S,s,H,h0,h,D,tini)
# 
# --- Périodisation ---
# 
N = 2
theta_periodise = np.tile(theta_T,N)   
w_periodise = np.append(w_T, T/tau + w_T)
# 
theta_pt_periodise = np.tile(theta_T_pt,N)
# 
nf = nf+1
plt.figure(nf)
plt.subplot(2,2,1)
plt.title('Oscillations de relaxation de la hauteur de la surface libre \n $s=%0.1fcm^2$'%s+', $D=%0.3f$L/s'%D)
plt.plot(w_periodise,theta_periodise)
plt.xlim([0,N*T/tau])
plt.grid(True,'major','both')
plt.xlabel('$t/\\tau$')
plt.ylabel('$z/h$')
plt.legend((['$\\epsilon=$%0.2f'%eps+ ', ' + '$D/D_0$=%0.2f'%r+ ', '+ '$D_1/D_0$=%0.2f'%etalim]))# 
plt.subplot(2,2,2)
plt.title("Portrait de phase de l'oscillation")
plt.plot(theta_periodise,theta_pt_periodise)
plt.xlabel('$\\theta = z/h$')
plt.ylabel('$\\mathrm{d}_w \\theta, \quad w = t/\\tau$')
plt.grid(True,'major','both')
# 
# #############################################################################
# 
S = 80#cm**2
s = 1#cm**2
# -----
H = 1/2#m
h0 = 3/5*H#m
h = 1.2#m
# 
D = 0.40#LL/s
# 
tini = 0
# 
# #############################################################################
# 
# --- Oscillation de relaxation
w_T, theta_T, theta_T_pt, T, tr, etalim, r, eps, tau = TantalOscillator(nf,Np,S,s,H,h0,h,D,tini)
# 
# --- Périodisation ---
# 
N = 2
theta_periodise = np.tile(theta_T,N)   
w_periodise = np.append(w_T, T/tau + w_T)
# 
theta_pt_periodise = np.tile(theta_T_pt,N)
# 
plt.subplot(2,2,3)
plt.plot(w_periodise,theta_periodise)
plt.title('Oscillations de relaxation de la hauteur de la surface libre \n $s=%0.1fcm^2$'%s+', $D=%0.3f$L/s'%D)
plt.xlim([0,N*T/tau])
plt.grid(True,'major','both')
plt.xlabel('$t/\\tau$')
plt.ylabel('$z/h$')
plt.legend((['$\\epsilon=$%0.2f'%eps+ ', ' + '$D/D_0$=%0.2f'%r+ ', '+ '$D_1/D_0$=%0.2f'%etalim]))
# 
plt.subplot(2,2,4)
plt.title("Portrait de phase de l'oscillation B")
plt.plot(theta_periodise,theta_pt_periodise)
plt.xlabel('$\\theta = z/h$')
plt.ylabel('$\\mathrm{d}_w \\theta, \quad w = t/\\tau$')
plt.grid(True,'major','both')