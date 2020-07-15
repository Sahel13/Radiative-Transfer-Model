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

data = np.genfromtxt('absCoef.csv', delimiter=',')
concentrations = [200, 420, 600]

plt.figure(figsize=(13, 9))
for ppmv in concentrations:
    wavelengths = []
    emissionTemperatures = []
    mixingRatio = CO2MixingRatio(ppmv)
    rhoExp = rho_0 * mixingRatio * np.exp(-z/H)

    for row in range(len(data)):
        absCoef = data[row][0]
        wavelength = data[row][1]

        B = np.ones(numPoints)
        for i in range(numPoints):
            B[i] = BBintensity(z[i], earthTempProfile, wavelength)

        def dtau(y, x):
            return -absCoef * earthDensityProfile(x) * mixingRatio
        
        dictInput = {
            'sunOde': dtau,
            'sunIntensity': 0,
            'zValues': z,
            'earthOde': dtau, # not necessary
            'earthIntensity': 0 # not necessary
        }

        profiler = vertProfile(dictInput)
        tau, _ = profiler.sunProfile(1)

        Q =  B * rhoExp * absCoef * np.exp(-tau)

        zMax = z[0]
        QMax = Q[0]
        for index in range(numPoints):
            if QMax < Q[index]:
                QMax = Q[index]
                zMax = z[index]

        emissionTemp = earthTempProfile(zMax)

        wavelength = wavelength * 1000000
        wavelengths.append(wavelength)
        emissionTemperatures.append(emissionTemp)

    plt.plot(wavelengths, emissionTemperatures, '-o', label=f"For {ppmv}ppmv")

plt.xlabel("Wavelength in micrometers", fontsize=13)
plt.ylabel("Temperature in K", fontsize=13)
plt.legend(fontsize=13)
plt.show()