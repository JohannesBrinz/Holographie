''' Versuch HO
 Johannes Brinz & Caterina Vanelli
 Datum: 10.01.2021
 '''

import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from scipy import optimize
import math
import matplotlib.font_manager as fm
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.signal import find_peaks


#Datenimmport
kohärenz = pd.DataFrame()

kohärenz["distance"] = [0.093, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
for i in range(0, 19):
    kohärenz["distance"][i+1] = kohärenz["distance"][i] + 0.025
    kohärenz["distance"][i] = kohärenz["distance"][i] - 0.093

kohärenz["intensity"] = [1, 0.9, 0.9, 0.65, 0.4, 0.15, 0.1, 0.05, 0, 0.05, 0.1, 0.2, \
                        0.3, 0.5, 0.7, 0.75, 0.9, 1, 1]

#Kontrastfunktion
tau = 0
Omega = 0
X = 0
def kontrast(x, N, f):
    tau = x / (3e8)
    Omega = (N-1)*2*np.pi * f
    return 1/N * abs(  np.sin( (N* Omega * tau)/(2*(N-1)) )  / np.sin( (Omega * tau)/(2*(N-1)) ) )

def kontrastII(x, N, f):
    tau = x / (3e8)
    X = tau * 2*np.pi * f
    return 1/6 + 1/3 * np.sqrt( 9/4 + 2*np.cos(X)*np.cos(2*X) + np.cos(X)*np.cos(3*X)\
        + np.cos(2*X)*np.cos(3*X) + 2*np.sin(X)*np.sin(2*X) + np.sin(X)*np.sin(3*X)\
            + np.sin(2*X)*np.sin(3*X))

#Eingabe
x_err = 0.001
I_err = 0.1

N = 4           #float(input("Moden-Zahl: "))               #Moden Zahl
f = 320e6        #float(input("Frequenzbreite: "))           #Frequenzbreite
t_c = 1/(N * f)

#fit

params, params_cov = optimize.curve_fit(kontrastII, kohärenz["distance"], kohärenz["intensity"], p0 = [4, 320e6])


print(t_c)
print(t_c * 3e8)
print(params)


#plot

plt.errorbar(kohärenz["distance"], kohärenz["intensity"], xerr = x_err, yerr = I_err,\
            linewidth = 2, fmt='+', color = "green", capsize=3)
plt.errorbar(np.linspace(0, 0.5, 1000), kontrastII(np.linspace(0,0.5, 1000), N, f))
plt.errorbar(np.linspace(0, 0.5, 1000), kontrast(np.linspace(0, 0.5, 1000), N, f), lw=1, fmt = "--", c = "black")
#plt.errorbar(np.linspace(0, 0.5, 1000), np.linspace((1/2.718),(1/2.718), 1000))
plt.title('Contrast function', fontsize = 15)
#plt.text( 0, 0,"$L_c = $" + str(round(t_c * 3*10e8, 3)))
plt.xlabel('$\Delta x$ [m]' , fontsize = 13)
plt.ylabel('K [-]', fontsize = 13)
plt.grid(True)
plt.legend(['experimental data', "weighted modes", "identical modes"], fontsize = 13)
plt.savefig('Plots/Kohärenzlänge.png', dpi=300)
plt.clf()


plt.errorbar(kohärenz["distance"], kohärenz["intensity"], xerr = x_err, yerr = I_err,\
            linewidth = 2, fmt='+', color = "green", capsize=3)
plt.errorbar(np.linspace(0, 0.5, 1000), kontrastII(np.linspace(0,0.5, 1000), params[0], params[1]))
plt.errorbar(np.linspace(0, 0.5, 1000), kontrast(np.linspace(0, 0.5, 1000), params[0], params[1]), lw=1, fmt = "--", c = "black")
plt.title('Fitted contrast function', fontsize = 15)
plt.text( 0.13, 0.65,"$\Delta\\nu_{mod} = $" + str(round(params[1]/1e6, 1)) + "MHz", fontsize = 13)
plt.xlabel('$\Delta x$ [m]' , fontsize = 13)
plt.ylabel('K [-]', fontsize = 13)
plt.grid(True)
plt.legend(['experimental data', "weighted modes", "identical modes"], fontsize = 13)
plt.savefig('Plots/Kohärenzlänge_fit.png', dpi=300)
plt.clf()

N=3

plt.errorbar(kohärenz["distance"], kohärenz["intensity"], xerr = x_err, yerr = I_err,\
            linewidth = 2, fmt='+', color = "green", capsize=3)
#plt.errorbar(np.linspace(0, 0.5, 1000), kontrastII(np.linspace(0,0.5, 1000), N, f))
plt.errorbar(np.linspace(0, 0.5, 1000), kontrast(np.linspace(0, 0.5, 1000), N, f), lw=1, fmt = "--", c = "black")
#plt.errorbar(np.linspace(0, 0.5, 1000), np.linspace((1/2.718),(1/2.718), 1000))
plt.title('Contrast function', fontsize = 15)
#plt.text( 0, 0,"$L_c = $" + str(round(t_c * 3*10e8, 3)))
plt.xlabel('$\Delta x$ [m]' , fontsize = 13)
plt.ylabel('K [-]', fontsize = 13)
plt.grid(True)
plt.legend(['experimental data', "identical modes"], fontsize = 13)
plt.savefig('Plots/Kohärenzlänge_3Moden.png', dpi=300)
plt.clf()
