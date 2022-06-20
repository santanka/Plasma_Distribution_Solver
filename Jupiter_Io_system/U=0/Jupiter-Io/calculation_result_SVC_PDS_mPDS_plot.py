import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
import numpy as np

data = np.genfromtxt(r'/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/U=0/Jupiter-Io/Jupiter_Io_calculation_result_SVC_PDS_mPDS.csv', delimiter=',', unpack=True)

B_ratio = data[0, :]
mlat_degree = data[1, :]
number_density_SVC = data[2, :]
number_density_PDS = data[3, :]
number_density_mPDS = data[4, :]
particle_flux_density_SVC = data[5, :]
particle_flux_density_PDS = data[6, :]
particle_flux_density_mPDS = data[7, :]
pressure_perp_SVC = data[8, :]
pressure_perp_PDS = data[9, :]
pressure_perp_mPDS = data[10, :]
pressure_para_SVC = data[11, :]
pressure_para_PDS = data[12, :]
pressure_para_mPDS = data[13, :]
pressure_dynamic_SVC = data[14, :]
pressure_dynamic_PDS = data[15, :]
pressure_dynamic_mPDS = data[16, :]
pressure_total_SVC = data[17, :]
pressure_total_PDS = data[18, :]
pressure_total_mPDS = data[19, :]

channel = 6

if (channel == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{n_{i}}{Z}$', yscale='log')
    ax.plot(mlat_degree, number_density_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, number_density_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, number_density_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()

if (channel == 2):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{n_{i} V_{\parallel i}}{Z v_{th}}$', yscale='log')
    ax.plot(mlat_degree, particle_flux_density_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, particle_flux_density_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, particle_flux_density_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()

if (channel == 3):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{P_{\perp i}}{Z T_{b}}$', yscale='log')
    ax.plot(mlat_degree, pressure_perp_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_perp_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_perp_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()

if (channel == 4):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{P_{\parallel i}}{Z T_{b}}$', yscale='log')
    ax.plot(mlat_degree, pressure_para_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_para_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_para_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()

if (channel == 5):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{P_{\mathrm{dynamic} \, i}}{Z T_{b}}$', yscale='log')
    ax.plot(mlat_degree, pressure_dynamic_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_dynamic_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_dynamic_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()

if (channel == 6):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{P_{\mathrm{total} \, i}}{Z T_{b}}$', yscale='log')
    ax.plot(mlat_degree, pressure_total_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_total_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, pressure_total_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()