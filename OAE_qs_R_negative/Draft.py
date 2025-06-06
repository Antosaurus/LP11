# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 22:45:09 2022
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
import scipy.interpolate as spinterp
import scipy.special as spsp
from scipy.optimize import fsolve
#
#-Plot figures
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.close('all')
mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.sans-serif'] = 'Arial'
#mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['figure.dpi'] = 110
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.grid.axis'] = 'both'
mpl.rcParams['axes.grid.which'] = 'minor'
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['figure.subplot.hspace']=  0.6
mpl.rcParams['lines.linewidth']  = 1.4
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['figure.figsize'] = 1.5*7.4, 1.5*5.8
# mpl.rcParams['text.usetex'] = True
from mpl_toolkits.mplot3d import Axes3D
#
"""
Début
"""
nf = 0
# 
"""
Constantes et paramètres
"""
# 
R2 = 1000
Rprime = 10*1000
R1 = 1000
Rn = R1*Rprime/R2
R = 100*1000
RnR = Rn/R
Vsat = 12
v0 = (R1/(R1+R2))*Vsat
# 
L = 100*0.001
C = 0.1*10**-6
# 
m = 0.5*(1/R)*np.sqrt(L/C)
print('m = %.3f'%m)
# 
omega0 = np.sqrt(1/(L*C))
print('omega0 = %.3f rad/s'%omega0)
# 
# 
"""
Calcul de l'amplitude du cycle limite avec une méthode grossiere
"""
# -- Approx de la caractéristique en N
# x=np.linspace(0,np.pi,1000)
# nf = nf+1
# plt.figure(nf)
# plt.plot(x,x*(x**2-1))
# 
nf = nf +1
plt.figure(nf)
x=np.linspace(-1,1,1000)
plt.plot(x,3*x**2-1)
plt.grid('True','major','both')
plt.title('Approximation de $\phi$ (dérivée de la caractéristique de D)')
plt.xlabel('$x/x_s$')
plt.ylabel('$\\phi_{app}$')
plt.xlim([-1,1])
plt.ylim([-1,2])
# 
u0_approx = Vsat*np.sqrt((4/3)*(1-RnR))
print('u_0_approx = %.3f V'%u0_approx)
# 
"""
Calcul exact de l'amplitude du cycle limite
"""
# 
A = 0.5 + (np.pi/4)*(1-RnR)
# 
a = (np.cos(A))**-1
u0 = v0*a
print('u_0 = %.3f V'%u0)
# 
