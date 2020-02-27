#!/usr/bin/python
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

qdots = [0.,.010,.025,.050,.075,.100,.150,.200,.250] # W

qdot_bg = 0.15 # W

down = False

capthex = 1.0 # K
rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8.5/1000. # m^3 -- volume of He-II bottle
m = 3. # dimensionless -- the GM exponent

dx = 0.01 # m -- "length" of channel (hole)
dy = dx # m -- side length of hole area
ahole = dy*dy # m^2 -- effective area of hole

def func(t,a,b,c):
    return a*np.exp(-t/b)+c

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

# given T_{HEX}, and geometry factors, calculate T_{bath} in equilibrium
def capt_infty(qdot):
    ftinv_integral=0.
    ftinv_goal=dx*(qdot/ahole)**m
    captprime=capthex
    deltacaptstep=0.0001
    while (ftinv_integral<ftinv_goal):
        ftinv_integral=ftinv_integral+deltacaptstep*ftinv_sciver(captprime)
        captprime=captprime+deltacaptstep
    return captprime
        
dt = .01 # s time step
nsteps = 100000

iq=0
captend=np.empty(len(qdots))
captendthy=np.empty(len(qdots))
fig1=plt.figure()
for qdot in qdots:

    # First, calculate starting temperature by assuming it has equilibrated.
    if(down):
        captstart=capt_infty(qdot_bg+qdot)
    else:
        captstart=capt_infty(qdot_bg)
    print('Starting temperature is %f'%captstart)

    t = 0 # time
    capt = captstart # temperature

    ts = np.empty(nsteps)
    capts = np.empty(nsteps)

    ftinv_integral=dx*((qdot_bg)/ahole)**m
    for i in np.arange(nsteps):
        ts[i]=t
        capts[i]=capt
        dcapt=dt*(qdot+qdot_bg-ahole*(ftinv_integral/dx)**(1/m))/(rho*v*c(capt))
        t=t+dt
        capt=capt+dcapt
        ftinv_integral=ftinv_integral+dcapt*ftinv_sciver(capt)
        #print(i,t,dcapt,capt,ftinv_integral)

    plt.plot(ts,capts,label='%4.3f W'%qdot)

    # fit curve
    popt,pcov=curve_fit(func,ts,capts)
    print(qdot)
    print(popt)
    #plt.plot(ts,func(ts,*popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    
    captend[iq]=capt

    # Calculate theoretical balance of qdot and ftinv integral

    ftinv_integral=dx*((qdot_bg)/ahole)**m
    ftinv_goal=dx*((qdot+qdot_bg)/ahole)**m
    captprime=captstart
    deltacaptstep=0.0001
    while (ftinv_integral<ftinv_goal):
        ftinv_integral=ftinv_integral+deltacaptstep*ftinv_sciver(captprime)
        captprime=captprime+deltacaptstep

    captendthy[iq]=captprime
    iq=iq+1

plt.xlabel('Time (s)')
plt.ylabel('Bottle temperature (K)')
plt.legend()

#for qdot in qdots:
fig2=plt.figure()
plt.plot(qdots,captend,label='after 1000 s')
plt.plot(qdots,captendthy,label='at infinity')
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Power (W)')
plt.ylabel('Final Temperature (K)')
plt.legend()
plt.show()
