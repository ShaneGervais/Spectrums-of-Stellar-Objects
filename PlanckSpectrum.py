"""
The Planck function describes the energy spectrum
of a black body with temperature T. It solves the 
energy emitted per unit surface area.

We plot the energies emitted for given stars
"""

import numpy as np
from scipy.constants import h, c, k
import astropy.units as units
from astropy.constants import R_sun, M_sun, L_sun
from math import pi
from scipy.constants import sigma #Stefan-Boltzmann constant
import matplotlib.pyplot as plt

# made by petrklus at https://gist.github.com/petrklus/b1f427accdf7438606a6
#used to get different colors for different temperatures
from rgb_to_kelvin import convert_K_to_RGB

#Calculates the luminosity of a star with the Stefan-Boltzmann law
def luminosity(R, Teff):
    
    #area
    A = 4*pi*R**2

    return A*sigma*Teff**4

#Calculate the radius of a star
def radius(L, Teff):
    return np.sqrt(L/(4*pi*sigma*Teff**4))

#arguments for each star
def stellar_parameters(*args):
    return {"R"     : args[0].to(units.m),
            "Teff"  : args[1].to(units.K),
            "M"     : args[2].to(units.kg),
            "L"     : args[3].to(units.W)}

#Stars we will use
stars = {
    'Bernard\'s Star':
        stellar_parameters(0.196*R_sun, 3.13e3*units.K, 0.144*M_sun, luminosity(0.196*R_sun, 3.13e3*units.K)*units.W/(units.K**4 *units.m**2)),
    'Sirius A':
        stellar_parameters(1.711*R_sun, 9.94e3*units.K, 2.06*M_sun, luminosity(1.711*R_sun, 9.94e3*units.K)*units.W/(units.K**4 *units.m**2)),
    'Sirius B':
        stellar_parameters(5.8e3*units.km, 2.48e4*units.K, 1.02*M_sun, luminosity(5.8e3*units.km, 2.48e4*units.K)*units.W/(units.K**4 *units.m**2)),
    'Arcturus':
        stellar_parameters(25.4*R_sun, 4.29e3*units.K, 1.1*M_sun, luminosity(25.4*R_sun, 4.29e3*units.K)*units.W/(units.K**4 *units.m**2)),
    'Betelgeuse':
        stellar_parameters(6.4e8*units.km, 3.59e3*units.K, 1.2*M_sun, luminosity(6.4e8*units.km, 3.59e3*units.K)*units.W/(units.K**4 *units.m**2)),
    'Aldebaran':
        stellar_parameters(radius(4.4e2*L_sun, 3.9e3*units.K)*((units.K**2 * units.m)/units.W**(1/2)), 3.9e3*units.K, 1.2*M_sun, 4.4e2*L_sun),
    'Bellatrix':
        stellar_parameters(radius(9.21e2*L_sun, 2.2e4*units.K)*((units.K**2 * units.m)/units.W**(1/2)), 2.2e4*units.K, 8.6*M_sun, 9.21e2*L_sun)
}

def planck_spectrum(wavelength, T):
    return 2*h*c**2/(wavelength**5 * (np.exp(h*c/(wavelength*k*T)) - 1))

#initialize array of temp
temperatures = np.zeros(len(stars) + 1) #this one will be sorted
effective_temp = np.zeros(len(stars) + 1) #this one will not be
lumins = np.zeros(len(stars) + 1)
radius = np.zeros(len(stars) + 1)
mass = np.zeros(len(stars) + 1)

#enum values
for i, key in enumerate(stars):
    temperatures[i] = stars[key]['Teff'].value
    effective_temp[i] = stars[key]['Teff'].value
    lumins[i] = stars[key]['L'].value
    radius[i] = stars[key]['R'].value
    mass[i] = stars[key]['M'].value


#initial values are the Sun
temperatures[-1] = 5778 #K
effective_temp[-1] = 5778 #K
lumins[-1] = L_sun/units.W
radius[-1] = R_sun/units.m
mass[-1] = M_sun/units.kg

#sorted temperatures
temperatures = np.sort(temperatures) 

#grids
n = 1000
#max wavelength
lambda_max = 2e-6
wavelength = np.linspace(lambda_max/n, lambda_max, n)

#plotting
plt.figure("Planck Spectrum")

for T in temperatures:
    color = tuple([val/255 for val in convert_K_to_RGB(T)])

    plt.semilogy(1e9*wavelength, 1e-12*planck_spectrum(wavelength, T), color=color, label="{:.0f} K".format(T))

plt.xlabel("$\lambda$ [nm]")
plt.xlim(0, 1e9*lambda_max)
plt.ylabel("$B_\lambda(T) $" + "[$\mathrm{kW\,m^{-2}\,nm^{-1}\, sr^{-1}}$]")
plt.ylim(0.1,5e4)
plt.legend(loc="upper right")
plt.show()


#uncomment this to analyze the data in terminal. May also be transfered to csv file
"""
print("T                M             L               R   ")
print("===================================================")

for i in range(len(stars)+1):
    print("%e %e %e %e\n" % (effective_temp[i], mass[i], lumins[i], radius[i]))
"""

plt.figure("Hertzsprung-Russell diagram")
plt.ylabel("Luminosity (/W)")
plt.xlabel("Temperature (/K)")
plt.gca().invert_xaxis()
for i in range(len(stars)):
    plt.loglog(effective_temp[i], lumins[i], 'o', color=tuple([val/255 for val in convert_K_to_RGB(effective_temp[i])]))
plt.show()