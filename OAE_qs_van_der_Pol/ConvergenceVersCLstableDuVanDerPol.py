# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:53:49 2022
Domaine : Mécanique 
Sous-Domaine : Oscillateurs ; portraits de phase et non linéarités
Chapitre : van der Pol
Appellation : Convergence vers le cycle limite stable du van der Pol
@author: Antonin
"""
#-Clean Workspace
get_ipython().magic('reset -sf')
#
#-Import libraries
import numpy as np
import scipy.interpolate as spinterp
import scipy.special as spsp
from scipy.integrate import solve_ivp
#
#-Plot figures
import matplotlib as mpl
import matplotlib.pyplot as plt
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
mpl.rcParams['lines.linewidth']  = 1.4
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['figure.figsize'] = 1.5*5.8, 1.5*5.8
# mpl.rcParams['text.usetex'] = True
from mpl_toolkits.mplot3d import Axes3D
#
"""
Début
"""
nf = 0
#
"""
Résolution numérique de l'équation de van der Pol :
    d_t**2{x} + eps*(x**2-1)*d_t{x} + x = 0
"""
# Nécessite "from scipy.integrate import solve_ivp".
# Il faut récrire l'équation sous la forme d'un système du 1er ordre.
# Ici : d_t{V} = f(V,t).
# Avec : V=[x,d_t{x}] et f:[a,b]->[b,-eps*(a**2-1)*b].
# Conditions initiales : V1 = [x1,0], x1 = a1*cos(theta1)
# 
# --- Parametre epsilon ---
# 
# 
# 
# --- Définition de la fonction f ---
# 
def f(t,V):
    return [V[1],-eps*(V[0]**2-1)*V[1]-V[0]]
# 
# --- Intervalle de temps ---
# 
eps = 0.10
Ncycles = 7*10
tstart = 0
Np = 10000
tend = Ncycles*np.pi
# 
t_span = [tstart,tend]
t = np.linspace(tstart, tend, Np)
# 
# --- Conditions initiales ---
# - cycle avec gain d'énergie
a1 = 0.15 #cycle limite à a = 2
theta1 = 0
x1 = a1*np.cos(theta1)
y1 = a1*np.sin(theta1)
# 
Vini = [x1, y1]
# 
# - cycle avec perte d'énergie
a2 = 3
theta2 = 0
x2 = a2*np.cos(theta2)
y2 = a2*np.sin(theta2)
# 
Vini2 = [x2,y2]
# 
# --- Résolution numérique ---
# 
V = solve_ivp(f, t_span, Vini, method='RK45', t_eval=t, dense_output=False, events=None, vectorized=False, args=None)
#
V2 = solve_ivp(f, t_span, Vini2, method='RK45', t_eval=t, dense_output=False, events=None, vectorized=False, args=None)
# 
#--- Cycle limite ---
# 
xlimth = 2*np.cos(t)
ylimth = -2*np.sin(t) 
# 
# --- Plot du portrait de phase ----
nf = nf+1
plt.figure(nf)
plt.plot(Vini[0],Vini[1],'o')
plt.plot(V.y[0],V.y[1],color='black')
plt.plot(xlimth,ylimth,color='red')
plt.plot(Vini2[0],Vini2[1],'o')
plt.plot(V2.y[0],V2.y[1],color='blue')
plt.title("Portrait de phase de l'oscillateur de van der Pol")
plt.grid(True,'major','both')
plt.xlabel('$x(t)$')
plt.ylabel('$y(t)=\\mathrm{d}_tx$')
plt.legend(['$(x(0),y(0))=(%0.2f,%0.2f)$'%(x1,y1),('$\\epsilon=%0.2f,$'%eps+' $a(0)=%0.2f$'%a1),'$\\epsilon=0$, $a(t)=a_\infty=2$','$(x(0),y(0))=(%0.2f,%0.2f)$'%(x2,y2),('$\\epsilon=%0.01f,$'%eps+' $a(0)=%0.2f$'%a2),])
plt.tight_layout()
# 
# --- Plot de la solution versus le temps
# 
# --- amplitude approchée (méthode Krylov et Bogolioubov)---
# 
aKB = 2/(np.sqrt(1-(1-4/a1**2)*np.exp(-eps*t)))
aKB2 = 2/(np.sqrt(1-(1-4/a2**2)*np.exp(-eps*t)))
# 
nf = nf+1
plt.figure(nf)
# 
plt.subplot(1,2,1)
plt.plot(t,V.y[0])
plt.plot(t,aKB,'--',color='purple')
plt.title(("Graphe de $x(t)$ (sol. num. exacte) et de $a_{KB}(t)$ calculée avec \n l'approximation de Krylov et Bogolioubov.\n $\\epsilon=%0.2f$"%eps+", $a(0)=$%0.2f"%a1))
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),a_{KB}(t)$')
plt.legend(['$x(t)$','$a_{KB}(t)$'])
# 
plt.subplot(1,2,2)
plt.title(("Graphe de $x(t)$ (sol. num. exacte) et de $a_{KB}(t)$ calculée avec \n l'approximation de Krylov et Bogolioubov.\n $\\epsilon=%0.2f$"%eps + ", $a(0)=$%0.2f"%a2))
plt.plot(t,V2.y[0])
# plt.plot(t,2*np.cos(t))
plt.plot(t,aKB2,'--',color='purple')
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),a_{KB}(t)$')
plt.legend(['$x(t)$','$a_{KB}(t)$'])
# 
plt.tight_layout()
# 
# --- Solution approchée analytique versus solution exacte numérique
# - thetaKB
thetaKB = theta1 - t
thetaKB2 = theta2 - t
# - xKB
xKB = aKB*np.cos(thetaKB)
xKB2 = aKB2*np.cos(thetaKB2)
# -xL
u = 0.5*eps*(1-a1**2/4)
v = np.sqrt(1-u**2)
xL = (a1/v)*np.exp(u*t)*(v*np.cos(v*t)-u*np.cos(v*t))
# 
u = 0.5*eps*(1-a2**2/4)
v = np.sqrt(1-u**2)
xL2 = (a2/v)*np.exp(u*t)*(v*np.cos(v*t)-u*np.cos(v*t))
# - Plot
nf = nf+1
plt.figure(nf)
# 
plt.subplot(2,1,1)
plt.plot(t,V.y[0])
plt.plot(t,xKB,'--',color='purple')
plt.title(("Graphe de $x(t)$ et de $x_{KB}(t)$ calculé avec \n l'approximation de Krylov et Bogolioubov.\n $\\epsilon=%0.2f$"%eps + ", $a(0)=$%0.2f"%a2+", $\\theta(0)=$%0.2f"%theta1))
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),x_{KB}(t)$')
plt.legend(['$x(t)$','$x_{KB}(t)$'])
# 
plt.subplot(2,1,2)
plt.plot(t,V2.y[0])
plt.plot(t,xKB2,'--',color='purple')
plt.title(("Graphe de $x(t)$ et de $x_{KB}(t)$ calculé avec \n l'approximation de Krylov et Bogolioubov.\n $\\epsilon=%0.2f$"%eps + ", $a(0)=$%0.2f"%a2+", $\\theta(0)=$%0.2f"%theta2))
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),x_{KB}(t)$')
plt.legend(['$x(t)$','$x_{KB}(t)$'])
# 
nf = nf+1
plt.figure(nf)
# 
plt.subplot(2,1,1)
plt.plot(t,V.y[0])
plt.plot(t,xL,'--',color='purple')
plt.title(("Graphe de $x(t)$ et de $x_{L}(t)$ calculé avec \n la pseudo-linéarisation de l'équation NL.\n $\\epsilon=%0.2f$"%eps + ", $a(0)=$%0.2f"%a2+", $\\theta(0)=$%0.2f"%theta1))
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),x_{L}(t)$')
plt.legend(['$x(t)$','$x_{L}(t)$'])
# 
plt.subplot(2,1,2)
plt.plot(t,V2.y[0])
plt.plot(t,xL2,'--',color='purple')
plt.title(("Graphe de $x(t)$ et de $x_{L}(t)$ calculé avec \n la pseudo-linéarisation de l'équation NL.\n $\\epsilon=%0.2f$"%eps + ", $a(0)=$%0.2f"%a2+", $\\theta(0)=$%0.2f"%theta2))
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t),x_{L}(t)$')
plt.legend(['$x(t)$','$x_{L}(t)$'])
# 
# --- Spectre du régime associé au cycle limite ---
# --- Extraction du régime établi
Tlimth = 2*np.pi
freqlimth = 1/Tlimth
# 
nf = nf+1
plt.figure(nf)
# - Cycle avec gain d'énergie
threshold = 1.99
listind = np.where(V.y[0]>=threshold)
indstartcl = np.min(listind)
xlimnum = V.y[0][indstartcl:-1]
vect = t[indstartcl:-1]
samplespacing = (t[-1]-t[indstartcl])/np.size(xlimnum)
Nptperiode = int(Tlimth / samplespacing)
# 
plt.subplot(2,2,1)
plt.plot(vect[-1-3*Nptperiode:-1],xlimnum[-1-3*Nptperiode:-1])
plt.plot(vect[-1-3*Nptperiode:-1],2*np.cos(vect[-1-3*Nptperiode:-1]),'--')
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t)$')
plt.legend(['$x(t)$','$x_0(t) = 2 \\mathrm{cos}(t)$'])
plt.title('Graphe de $x(t)$ en régime établi \n $\\epsilon=%0.2f$'%eps+' et $a(0)=%0.2f$'%a1)
# 
plt.subplot(2,2,2)
tfx = np.fft.rfft(xlimnum, norm='ortho')
freq = np.fft.rfftfreq(np.size(xlimnum),samplespacing)
plt.plot(freq,np.abs(tfx))
plt.axvline(x=freqlimth, color='brown')
plt.xlim([0,10*freqlimth])
plt.grid(True,'major','both')
plt.xlabel('$\\nu[Hz]$')
plt.ylabel('$\\left|\\tilde{x}(\\nu)\\right|$')
plt.legend(['FFT','$f_\infty=1/2\\pi$'])
plt.title('Module de la TF de $x(t)$ en régime établi')
# 
# - Cycle avec perte d'énergie
# 
# Attention si a2 est trop grand (on démarre trop loin du cycle limite)
#  un déphasage apparaît dans le régime établi par rapport au cosinus du régime 
#  epsilon = 0. 
# 
threshold = 2.01
listind = np.where(V2.y[0]>=threshold)
indstartcl = np.max(listind)
xlimnum = V2.y[0][indstartcl:-1]
vect = t[indstartcl:-1]
samplespacing = (t[-1]-t[indstartcl])/np.size(xlimnum)
Nptperiode = int(Tlimth / samplespacing)
#
plt.subplot(2,2,3)
plt.plot(vect[-1-3*Nptperiode:-1],xlimnum[-1-3*Nptperiode:-1])
plt.plot(vect[-1-3*Nptperiode:-1],2*np.cos(vect[-1-3*Nptperiode:-1]),'--')
plt.grid(True,'major','both')
plt.xlabel('$t$')
plt.ylabel('$x(t)$')
plt.legend(['$x(t)$','$x_0(t) = 2 \\mathrm{cos}(t)$'])
plt.title('Graphe de $x(t)$ en régime établi \n $\\epsilon=%0.2f$'%eps+' et $a(0)=%0.2f$'%a2)
# 
plt.subplot(2,2,4)
tfx = np.fft.rfft(xlimnum, norm='ortho')
Tlimth = 2*np.pi
freqlimth = 1/Tlimth
samplespacing = (t[-1]-t[indstartcl])/np.size(xlimnum)
freq = np.fft.rfftfreq(np.size(xlimnum),samplespacing)
plt.plot(freq,np.abs(tfx))
plt.axvline(x=freqlimth,color='brown')
plt.xlim([0,10*freqlimth])
plt.grid(True,'major','both')
plt.xlabel('$\\nu[Hz]$')
plt.ylabel('$\\left|\\tilde{x}(\\nu)\\right|$')
plt.legend(['FFT','$f_\infty=1/2\\pi$'])
plt.title('Module de la TF de $x(t)$ en régime établi')
# 

