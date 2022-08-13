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

channel = 8

if (channel == 1):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', yscale='log')
    ax.set_ylabel(r'$\frac{n_{i}}{Z}$', rotation=0, va='center', size=90)
    ax.plot(mlat_degree, number_density_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, number_density_PDS, label=r'PDS-1', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, number_density_mPDS, label=r'PDS-2', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=60)
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

if (channel == 7):
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.usetex'] = True
    plt.rcParams["font.size"] = 25

    fig = plt.figure(tight_layout=True)
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0], xlabel=r'MLAT [degree]', yscale='log')
    ax2 = fig.add_subplot(gs[1, 0], xlabel=r'MLAT [degree]', yscale='log')
    ax3 = fig.add_subplot(gs[2, 0], xlabel=r'MLAT [degree]', yscale='log')
    ax1.set_ylabel(r'$\frac{n_{i} V_{\parallel i}}{B_{i}} \frac{B_{b}}{Z v_{\mathrm{th}}}$', rotation=0, size=40)
    ax2.set_ylabel(r'$\frac{P_{\perp i}}{Z T}$', rotation=0, va='center', size=40)
    ax3.set_ylabel(r'$\frac{P_{\parallel i}}{Z T}$', rotation=0, size=40)

    ax1.set_title(r'(a)', x=-0.1, y=0.95)
    ax1.plot(mlat_degree, particle_flux_density_SVC  / B_ratio, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree, particle_flux_density_PDS / B_ratio, label=r'PDS-1', linestyle='solid', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree, particle_flux_density_mPDS  / B_ratio, label=r'PDS-2', linestyle='solid', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both", alpha=0.3)
    ax1.legend(fontsize=30)

    ax2.set_title(r'(b)', x=-0.1, y=0.95)
    ax2.plot(mlat_degree, pressure_perp_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree, pressure_perp_PDS, label=r'PDS-1', linestyle='solid', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree, pressure_perp_mPDS, label=r'PDS-2', linestyle='solid', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which='both', alpha=0.3)
    ax2.legend(fontsize=30)

    ax3.set_title(r'(c)', x=-0.1, y=0.95)
    ax3.plot(mlat_degree, pressure_para_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree, pressure_para_PDS, label=r'PDS-1', linestyle='solid', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree, pressure_para_mPDS, label=r'PDS-2', linestyle='solid', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which='both', alpha=0.3)
    ax3.legend(fontsize=30)

    plt.tight_layout()
    plt.show()

if (channel == 8):
    fig = plt.figure()
    plt.rcParams["font.size"] = 60
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']
    ax = fig.add_subplot(111, xlabel=r'MLAT [degree]', ylabel=r'$\frac{V_{\parallel i}}{v_{th}}$', yscale='log')
    ax.plot(mlat_degree, particle_flux_density_SVC / number_density_SVC, label=r'SVC', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, particle_flux_density_PDS / number_density_PDS, label=r'PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.plot(mlat_degree, particle_flux_density_mPDS / number_density_mPDS, label=r'm-PDS', linestyle='solid', linewidth='4', alpha=0.7)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.5)
    ax.legend(fontsize=40)
    plt.tight_layout()
    plt.show()
