import numpy as np
from numInt import *
import matplotlib.pyplot as plt

toa = 120 * 10**(3) # m
z = np.arange(0, toa, 1000)
numPoints = len(z)
zkm = z/1000
Cp = 1000 # J/kg
rho_0 = 1.225
H = 8500 # m
np.seterr('raise')

shortwave = 0.5 * 10**(-6)
longwave = 15 * 10**(-6)

klong = 10
kshort = 0.005

I_Sshort = BBsingleIntensity(5800, shortwave)
I_Slong = BBsingleIntensity(5800, longwave)

I_Eshort = BBsingleIntensity(290, shortwave)
I_Elong = BBsingleIntensity(290, longwave)

mixingRatio = CO2MixingRatio(200)
rhoExp = rho_0 * mixingRatio * np.exp(-z/H)

def sunOde(y, x):
    return (y - BBintensity(x, earthTempProfile, longwave)) * klong * earthDensityProfile(x) * mixingRatio

def earthOde(y, x):
    return (BBintensity(x, earthTempProfile, longwave) - y) * klong * earthDensityProfile(x) * mixingRatio

dictInput = {
            'sunOde': sunOde,
            'sunIntensity': I_Slong,
            'zValues': z,
            'earthOde': earthOde,
            'earthIntensity': I_Elong
        }

profiler = vertProfile(dictInput)
I, dIdz, sunI, sunDi, earthI, earthDi = profiler.twoWay(1)
dTdt = dIdz/(rhoExp*Cp)
dTdtmax = np.amax(np.abs(dTdt))

plt.figure()
plt.plot(dTdt/dTdtmax, zkm, label="For 200ppmv")

mixingRatio = CO2MixingRatio(420)
rhoExp = rho_0 * mixingRatio * np.exp(-z/H)

def sunOde(y, x):
    return (y - BBintensity(x, earthTempProfile, longwave)) * klong * earthDensityProfile(x) * mixingRatio

def earthOde(y, x):
    return (BBintensity(x, earthTempProfile, longwave) - y) * klong * earthDensityProfile(x) * mixingRatio

dictInput = {
            'sunOde': sunOde,
            'sunIntensity': I_Slong,
            'zValues': z,
            'earthOde': earthOde,
            'earthIntensity': I_Elong
        }

profiler = vertProfile(dictInput)
I, dIdz, sunI, sunDi, earthI, earthDi = profiler.twoWay(1)
dTdt = dIdz/(rhoExp*Cp)
dTdtmax = np.amax(np.abs(dTdt))

plt.plot(dTdt/dTdtmax, zkm, label="For 420ppmv")
plt.legend()
plt.xlabel("dT/dt (Normalized scale)")
plt.ylabel("Height in km")
plt.show()
