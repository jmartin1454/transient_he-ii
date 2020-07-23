#!/usr/bin/python
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from math import *

capthex = 1.0 # K
rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8.5/1000. # m^3 -- volume of He-II bottle
m = 3. # dimensionless -- the GM exponent

d = 0.15 # m -- diameter of He-II channel
A = pi*d**2/4 # m -- area of He-II channel

# van Sciver's fit for c(T)
def c(capt): # J/(kg*K)
    #return 100. # previous model
    if (capt<=.6):
        return 20.4*capt**3 # van Sciver Eq. (6.28)
    elif (capt<1.1):
        return 108.*capt**6.7 # van Sciver Eq. (6.29a)
    else:
        return 117.*capt**5.6 # van Sciver Eq. (6.29b)

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

nx=10 # number of position steps
L=1. # m, length of channel
tempi=1.0 # K, initial temperature
temp=[tempi]*nx # array of temperatures
dx=L/nx
x=[i*dx+dx/2 for i in range(nx)]

dt = .001 # s time step
nsteps = 200001
t = 0 # time

Qin=10.0 # W, heat flow in.
qin=Qin/A # W/m^2, heat flux in.

k=5000 # W/(m*K), for constant thermal conductivity model.


savetimes=[0.01,0.5,1.0,2.0,10.0,20.0]

for i in np.arange(nsteps):
    for j in np.arange(nx-1):
        if(j==0):
            temp[0]=temp[1]+qin*dx/k
            #print("Temp0 %f"%temp[0])
        else:
            leftderiv=k*(temp[j]-temp[j-1])/dx
            rightderiv=k*(temp[j+1]-temp[j])/dx
            secondderiv=(leftderiv-rightderiv)/dx
            dtemp=-dt/(rho*c(temp[j]))*secondderiv
            temp[j]=temp[j]+dtemp
            #print(j,leftderiv,rightderiv,dx,secondderiv,dtemp,temp[j])
    for savetime in savetimes:
        if(savetime-dt/2<t<savetime+dt/2):
            print(t,temp)
            plt.plot(x,temp,label='%4.2f s'%savetime)
    t=t+dt
plt.xlabel('Position (m)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

