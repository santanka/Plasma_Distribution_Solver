from cmath import sqrt
from struct import unpack
from tkinter import N
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
import numpy as np

#dipole magnetic field (L=5.85)
latitude_ionosphere = 65.1E0
grid_number = 10000     #integer

latitude_deg = np.linspace(0., latitude_ionosphere, grid_number)

latitude_rad = np.deg2rad(latitude_deg)

def magnetic_flux_density_ratio(latitude_rad):
    B_ratio = np.sqrt(1E0 + 3E0*np.sin(latitude_rad)**2E0) / np.cos(latitude_rad)**6E0
    return B_ratio

alpha = magnetic_flux_density_ratio(np.deg2rad(latitude_ionosphere))

beta = magnetic_flux_density_ratio(latitude_rad)

def number_density(alpha, beta):
    x = beta / alpha
    n = 5E-1 * (1E0 + np.sqrt(1E0 - x))
    return n

def mean_velocity(alpha, beta):
    x = beta / alpha
    Vpara = 1E0 / np.sqrt(np.pi) / (1E0 + np.sqrt(1E0 - x))
    return Vpara

def pressure_perpendicular(alpha, beta):
    x = beta / alpha
    Pperp = 2.5E-1 * (2E0 + (2E0 + x)*np.sqrt(1E0 - x))
    return Pperp

def pressure_parallel(alpha, beta):
    x = beta / alpha
    Ppara = 5E-1 * (1E0 + (1E0 + x)*np.sqrt(1E0 - x) + 2E0 / np.pi * (1E0 - 2E0*x) / (1E0 + np.sqrt(1E0 - x)))
    return Ppara

def temperature_perpendicular(number_density, pressure_perpendicular):
    Tperp = pressure_perpendicular / number_density
    return Tperp

def temperature_parallel(number_density, pressure_parallel):
    Tpara = pressure_parallel / number_density
    return Tpara

num = number_density(alpha, beta)
Vpara = mean_velocity(alpha, beta)
Pperp = pressure_perpendicular(alpha, beta)
Ppara = pressure_parallel(alpha, beta)
Tperp = temperature_perpendicular(num, Pperp)
Tpara = temperature_parallel(num, Ppara)

fig = plt.figure()
plt.rcParams["font.size"] = 60
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
ax = fig.add_subplot(111, xlabel=r'MLAT [degree]')
ax.plot(latitude_deg, num, label=r'$\frac{n_{i}}{Z}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.plot(latitude_deg, Vpara, label=r'$\frac{V_{\parallel i}}{v_{\mathrm{th}}}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.plot(latitude_deg, Pperp, label=r'$\frac{P_{\perp i}}{Z T}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.plot(latitude_deg, Ppara, label=r'$\frac{P_{\parallel i}}{Z T}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.plot(latitude_deg, Tperp, label=r'$\frac{T_{\perp i}}{T}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.plot(latitude_deg, Tpara, label=r'$\frac{T_{\parallel i}}{T}$', linestyle='solid', linewidth='4', alpha=0.7)
ax.minorticks_on()
ax.grid(which='both', alpha=0.5)
ax.legend(fontsize=40)
plt.tight_layout()
plt.show()

