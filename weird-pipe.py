#!/usr/bin/python
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from math import *

capthex = 1.0 # K
rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8.5/1000. # m^3 -- volume of He-II bottle
m = 3. # dimensionless -- the GM exponent

# van Sciver's fit for f(T)^{-1}
def ftinv_sciver(capt):
    t_lambda=2.17 # K
    smallt=capt/t_lambda # dimensionless
    multiplier=(smallt**5.7*(1.-smallt**5.7))**3 # dimensionless
    s_lambda=1559. # J/(kg K)
    A_lambda=1450. # (m s)/kg
    rho=145. # kg/m^3
    g_lambda=rho**2*s_lambda**4*t_lambda**3/A_lambda # W^3/(K m^5)
    return g_lambda*multiplier

capQ=10 # W
x=0
dx=0.001 # m
r=0.15/2 # m
A=pi*r**2

def A_changing(x):
    r_main=0.15/2
    r_small=0.13/2
    l_small=0.10
    c_small=1.0
    if (c_small-l_small/2<x<c_small+l_small/2):
        A=pi*r_small**2
    else:
        A=pi*r_main**2
    return A
        
numsteps=1500

x_save = np.empty(numsteps)
T_changing_save = np.empty(numsteps)
T_notchanging_save = np.empty(numsteps)

T_changing=1 # K
T_notchanging=1 # K

for i in range(0,numsteps):
    x=x+dx
    dT_changing=capQ**3/A_changing(x)**3*dx/ftinv_sciver(T_changing)
    dT_notchanging=capQ**3/A**3*dx/ftinv_sciver(T_notchanging)
    T_changing=T_changing+dT_changing
    T_notchanging=T_notchanging+dT_notchanging
    x_save[i]=x
    T_changing_save[i]=T_changing
    T_notchanging_save[i]=T_notchanging
    #print(x,T)

plt.plot(x_save,T_notchanging_save,label='notchanging')
plt.plot(x_save,T_changing_save,label='changing')
plt.legend()
plt.show()
