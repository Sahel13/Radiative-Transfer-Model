import numpy as np
import math

def BBintensity(height, tempProfile, wavelength):
	"""
	Function to return black body radiation intensity at a given height for a given temperature profile.
	Takes in arguments height, temperature profile function and wavelength.
	"""
	h = 6.62607015 * 10**(-34) # Planck's constant
	kb = 1.38064852 * 10**(-23) # Boltzman's constant
	c = 299792458 # Speed of light
	termTwo = math.exp(h*c/(wavelength*kb*tempProfile(height))) - 1
	termOne = 2*h*c**2/wavelength**5
	return termOne/termTwo

def BBsingleIntensity(temperature, wavelength):
	"""
	Returns blackbody radiation intensity given temperature and wavelength.
	"""
	h = 6.62607015 * 10**(-34) # Planck's constant
	kb = 1.38064852 * 10**(-23) # Boltzman's constant
	c = 299792458 # Speed of light
	termTwo = math.exp(h*c/(wavelength*kb*temperature)) - 1
	termOne = 2*h*c**2/wavelength**5
	return termOne/termTwo

def differentiate(func, x, h):
	if (x - h >= 0):
		return (func(x+h)-func(x-h))/(2*h)
	else:
		return (func(x+h)-func(x))/h

def earthTempProfile(x):
	"""
	Function that models earth's temperature profile.
	Returns temperature at altitude z with 0m <= z <= 120000m.
	"""
	if x <= 11000:
		return (290 - (70*x/11000))
	elif x <= 20000:
		return 220
	elif x <= 50000:
		return (220 + (50*(x-20000)/30000))
	elif x <= 90000:
		return (270 - (80*(x-50000)/40000))
	elif x <= 120000:
		return (190 + (170*(x-90000)/30000))
	else:
		return 0

def earthDensityProfile(x):
	"""
	Returns atmospheric density in kg/m3 at altitude z.
	"""
	return 1.225 * math.exp(-x/8500)

def CO2MixingRatio(ppmv):
	"""
	Returns mixing ratio of CO2 given its concentration in ppmv.
	"""
	return ppmv * 10**(-6) * 44.01/28.97

class vertProfile():
	"""
	Takes a dictionary as input with keys 'sunOde', 'earthOde', 'sunIntensity', 'earthIntensity' and 'zvalues'.
	"""
	def __init__(self, dictInput):
		self.sunOde = dictInput['sunOde']
		self.I_S = dictInput['sunIntensity']
		self.z = dictInput['zValues']
		self.earthOde = dictInput['earthOde']
		self.I_E = dictInput['earthIntensity']
		
	def intensityTopDown(self, low_limit, up_limit, boundary_value, step_size):
	    integral = boundary_value
	    x = up_limit
	    while (x >= low_limit):
	        integral -= self.sunOde(integral, x)*step_size
	        x -= step_size
	    return integral
	
	def intensityBottomUp(self, low_limit, up_limit, boundary_value, step_size):
	    integral = boundary_value
	    x = low_limit
	    while (x <= up_limit):
	        integral += self.earthOde(integral, x)*step_size
	        x += step_size
	    return integral

	def sunProfile(self, step_size):
		arrayLength = len(self.z)
		index = arrayLength - 1
		sunLoc = self.z[index]
		IArray = np.ones(arrayLength)
		dIdzArray = np.ones(arrayLength)
		IArray[index] = self.intensityTopDown(self.z[index], sunLoc, self.I_S, step_size)
		dIdzArray[index] = self.sunOde(IArray[index], self.z[index])
		while index > 0:
			index += -1
			intVal = self.intensityTopDown(self.z[index], self.z[index + 1], IArray[index + 1], step_size)
			IArray[index] = intVal
			dIdzArray[index] = self.sunOde(intVal, self.z[index])
		return (IArray, dIdzArray)

	def earthProfile(self, step_size):
		arrayLength = len(self.z)
		index = 0
		IArray = np.ones(arrayLength)
		dIdzArray = np.ones(arrayLength)
		IArray[index] = self.I_E
		dIdzArray[index] = self.earthOde(IArray[index], self.z[index])
		while index < arrayLength - 1:
			index += 1
			intVal = self.intensityBottomUp(self.z[index - 1], self.z[index], IArray[index - 1], step_size)
			IArray[index] = intVal
			dIdzArray[index] = self.earthOde(intVal, self.z[index])
		return (IArray, dIdzArray)
	
	def twoWay(self, step_size):
		sunI, sunDI = self.sunProfile(step_size)
		earthI, earthDI = self.earthProfile(step_size)
		I = np.subtract(sunI, earthI)
		DI = np.subtract(sunDI, earthDI)
		return (I, DI, sunI, sunDI, earthI, earthDI)