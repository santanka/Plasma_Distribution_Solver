import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_050_BC_004_min_074.csv", delimiter=',', unpack=True)

channel = 9
#1: number density, 2: plasma pressure perpendicular, 3: plasma pressure parallel, 4: dynamic pressure, 5: current density
#6: temperature perpendicular, 7: temperature parallel, 8: electrostatic potential, 9: pressure 3 kind, 
#10: compare to Knight relation (1973), 11: for paper, 12: Alfven wave transit time, 13: plasma beta 3 kind

series_number = 10

coordinate_FA = data[0, :]
length2planet = data[1, :]
mlat_rad = data[2, :]
mlat_degree = data[3, :]
magnetic_flux_density = data[4, :]
initial_electrostatic_potential = data[5, :]
electrostatic_potential = data[6, :]
number_density = data[7:series_number + 7, :]
charge_density = data[series_number + 7, :]
charge_density_Poisson = data[series_number + 8, :]
convergence_number = data[series_number + 9, :]
particle_flux_density = data[series_number + 10 : 2 * series_number + 10, :]
parallel_mean_velocity = data[2 * series_number + 10 : 3 * series_number + 10, :]
pressure_perp = data[3 * series_number + 10 : 4 * series_number + 10, :]
pressure_para = data[4 * series_number + 10 : 5 * series_number + 10, :]
pressure_dynamic = data[5 * series_number + 10 : 6 * series_number + 10, :]
temperature_perp = data[6 * series_number + 10 : 7 * series_number + 10, :]
temperature_para = data[7 * series_number + 10 : 8 * series_number + 10, :]
Alfven_speed = data[8 * series_number + 10, :]
Alfven_speed_per_lightspeed = data[8 * series_number + 11, :]
ion_inertial_length = data[8 * series_number + 12, :]
electron_inertial_length = data[8 * series_number + 13, :]
ion_Larmor_radius = data[8 * series_number + 14, :]
ion_acoustic_gyroradius = data[8 * series_number + 15, :]
electron_Larmor_radius = data[8 * series_number + 16, :]
current_density = data[8 * series_number + 17, :]

elementary_charge = 1.602176634E-19

if (channel == 1):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Number Density [$\mathrm{cm}^{-3}$]', yscale=r'log', ylim=(1E-2, 1E5))
    ax.plot(mlat_degree, number_density[0, :]*1E-6 + number_density[2, :]*1E-6, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, number_density[1, :]*1E-6 + number_density[3, :]*1E-6, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[4, :]*1E-6, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[5, :]*1E-6, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[6, :]*1E-6, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[7, :]*1E-6, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[8, :]*1E-6, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, number_density[9, :]*1E-6, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)


if (channel == 2):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Plasma Pressure (perpendicular) [$\mathrm{nPa}$]', yscale=r'log', ylim=(1E-5, 1E2))
    ax.plot(mlat_degree, (pressure_perp[0, :] + pressure_perp[2, :]) * 1E9, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, (pressure_perp[1, :] + pressure_perp[3, :]) * 1E9, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[4, :] * 1E9, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[5, :] * 1E9, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[6, :] * 1E9, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[7, :] * 1E9, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[8, :] * 1E9, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_perp[9, :] * 1E9, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)

if (channel == 3):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Plasma Pressure (parallel) [$\mathrm{nPa}$]', yscale=r'log', ylim=(1E-5, 1E2))
    ax.plot(mlat_degree, (pressure_para[0, :] + pressure_para[2, :]) * 1E9, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, (pressure_para[1, :] + pressure_para[3, :]) * 1E9, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[4, :] * 1E9, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[5, :] * 1E9, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[6, :] * 1E9, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[7, :] * 1E9, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[8, :] * 1E9, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_para[9, :] * 1E9, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)

if (channel == 4):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Dynamic Pressure (parallel) [$\mathrm{nPa}$]', yscale=r'log', ylim=(1E-5, 1E2))
    ax.plot(mlat_degree, (pressure_dynamic[0, :] + pressure_dynamic[2, :]) * 1E9, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, (pressure_dynamic[1, :] + pressure_dynamic[3, :]) * 1E9, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[4, :] * 1E9, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[5, :] * 1E9, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[6, :] * 1E9, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[7, :] * 1E9, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[8, :] * 1E9, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, pressure_dynamic[9, :] * 1E9, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)

if (channel == 5):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Current density [$\mathrm{\mu A / m^{2}}$]') #, yscale=r'log'
    ax.plot(mlat_degree, (current_density*1E6), c=r'blue', linestyle=r'solid', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)

if (channel == 6):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Temperature (perpendicular) [$\mathrm{eV}$]', yscale=r'log')
    ax.plot(mlat_degree, (pressure_perp[0, :] + pressure_perp[2, :]) / (number_density[0, :] + number_density[2, :]) / elementary_charge, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, (pressure_perp[1, :] + pressure_perp[3, :]) / (number_density[1, :] + number_density[3, :]) / elementary_charge, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[4, :] / elementary_charge, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[5, :] / elementary_charge, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[6, :] / elementary_charge, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[7, :] / elementary_charge, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[8, :] / elementary_charge, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_perp[9, :] / elementary_charge, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)

if (channel == 7):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'South ← MLAT [degree] → North', ylabel=r'Temperature (parallel) [$\mathrm{eV}$]', yscale=r'log')
    ax.plot(mlat_degree, (pressure_para[0, :] + pressure_para[2, :]) / (number_density[0, :] + number_density[2, :]) / elementary_charge, c=r'dimgrey', label=r'$\mathrm{H^{+}}$(Jupiter)', linestyle=r'solid', linewidth=r'4', alpha=0.7)
    ax.plot(mlat_degree, (pressure_para[1, :] + pressure_para[3, :]) / (number_density[1, :] + number_density[3, :]) / elementary_charge, c=r'blue', label=r'$\mathrm{e^{-}}$(Jupiter)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[4, :] / elementary_charge, c=r'purple', label=r'$\mathrm{H^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[5, :] / elementary_charge, c=r'orange', label=r'$\mathrm{O^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[6, :] / elementary_charge, c=r'green', label=r'$\mathrm{S^{+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[7, :] / elementary_charge, c=r'lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[8, :] / elementary_charge, c=r'deepskyblue', label=r'cold $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.plot(mlat_degree, temperature_para[9, :] / elementary_charge, c=r'hotpink', label=r'hot $\mathrm{e^{-}}$(Io)', linestyle=r'dotted', linewidth=r'4')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=r'upper left', borderaxespad=0, fontsize=40)

if (channel == 8):
    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2) = plt.subplots(2, 1) #, sharex=True
    fig.suptitle(r'Electrostatic Potential')
    ax1.set_xlabel(r'South ← MLAT [degree] → North')
    ax1.set_ylabel(r'[kV]')
    ax1.plot(mlat_degree, electrostatic_potential/1E3, linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both")
    ax2.set_xlabel(r'South ← MLAT [degree] → North')
    ax2.set_ylabel(r'[V]')
    ax2.plot(mlat_degree, electrostatic_potential, linewidth='4')
    ax2.set_ylim(-20, 5)
    ax2.minorticks_on()
    ax2.grid(which="both")
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 9):
    electron_perpendicular_pressure = pressure_perp[1, :] + pressure_perp[3, :] + pressure_perp[8, :] + pressure_perp[9, :]
    ion_perpendicular_pressure = pressure_perp[0, :] + pressure_perp[2, :] + pressure_perp[4, :] + pressure_perp[5, :] + pressure_perp[6, :] + pressure_perp[7, :]
    all_perpendicular_pressure = electron_perpendicular_pressure + ion_perpendicular_pressure

    electron_parallel_pressure = pressure_para[1, :] + pressure_para[3, :] + pressure_para[8, :] + pressure_para[9, :]
    ion_parallel_pressure = pressure_para[0, :] + pressure_para[2, :] + pressure_para[4, :] + pressure_para[5, :] + pressure_para[6, :] + pressure_para[7, :]
    all_parallel_pressure = electron_parallel_pressure + ion_parallel_pressure

    electron_total_pressure = electron_parallel_pressure / 3E0 + electron_perpendicular_pressure * 2E0 / 3E0
    ion_total_pressure = ion_parallel_pressure / 3E0 + ion_perpendicular_pressure * 2E0 / 3E0
    all_total_pressure = electron_total_pressure + ion_total_pressure

    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.set_title(r'perpendicular')
    ax1.set_xlabel(r'S ← MLAT [degree] → N')
    ax1.set_ylabel(r'[nPa]')
    ax1.set_yscale('log')
    ax1.set_ylim([1E-5, 1E2])
    ax1.plot(mlat_degree, electron_perpendicular_pressure*1E9, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree, ion_perpendicular_pressure*1E9, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree, all_perpendicular_pressure*1E9, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")

    ax2.set_title(r'parallel')
    ax2.set_xlabel(r'S ← MLAT [degree] → N')
    ax2.set_ylabel(r'[nPa]')
    ax2.set_yscale('log')
    ax2.set_ylim([1E-5, 1E2])
    ax2.plot(mlat_degree, electron_parallel_pressure*1E9, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree, ion_parallel_pressure*1E9, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree, all_parallel_pressure*1E9, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")

    ax3.set_title(r'total')
    ax3.set_xlabel(r'S ← MLAT [degree] → N')
    ax3.set_ylabel(r'[nPa]')
    ax3.set_yscale('log')
    ax3.set_ylim([1E-5, 1E2])
    ax3.plot(mlat_degree, electron_total_pressure*1E9, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree, ion_total_pressure*1E9, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree, all_total_pressure*1E9, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both")

    plt.subplots_adjust(wspace=0.4, hspace=0.6)

if (channel == 10):
    length = int((1+len(mlat_degree))/2)
    mlat_degree_half = mlat_degree[length-1:]

    electron_cold_temperature = (pressure_para[8, :] + pressure_perp[8, :]*2E0)/3E0 / number_density[8, :]
    electron_hot_temperature = (pressure_para[9, :] + pressure_perp[9, :]*2E0)/3E0 / number_density[9, :]
    electron_cold_temperature_boundary = electron_cold_temperature[length-1]
    electron_cold_temperature_end = electron_cold_temperature[len(mlat_degree)-1]
    electron_hot_temperature_boundary = electron_hot_temperature[length-1]
    electron_hot_temperature_end = electron_hot_temperature[len(mlat_degree)-1]

    electron_cold_number_density_end = number_density[8, len(mlat_degree)-1]
    electron_hot_number_density_end = number_density[9, len(mlat_degree)-1]
    electron_cold_number_density_boundary = number_density[8, length-1]
    electron_hot_number_density_boundary = number_density[9, length-1]

    magnetic_flux_density_half = magnetic_flux_density[length-1:]
    magnetic_flux_density_boundary = magnetic_flux_density[length-1]
    magnetic_flux_density_end = magnetic_flux_density[len(mlat_degree)-1]

    electrostatic_potential_boundary = electrostatic_potential[length-1]
    electrostatic_potential_end = electrostatic_potential[len(mlat_degree)-1]

    electron_mass = 9.1093837015E-031
    elementary_charge = 1.602176634E-19

    current_density_half = current_density[length-1:]

    current_density_Knight = electron_cold_number_density_boundary * np.sqrt(electron_cold_temperature_boundary / 2E0 / np.pi / electron_mass) \
        * (1E0 - (1E0 - magnetic_flux_density_boundary / magnetic_flux_density_end) \
            * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_cold_temperature_boundary) \
                * (magnetic_flux_density_end / magnetic_flux_density_boundary - 1E0))

    current_density_Knight = current_density_Knight \
        - electron_cold_number_density_end * np.sqrt(electron_cold_temperature_end / 2E0 / np.pi / electron_mass) \
            * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_cold_temperature_end) \
                * (1E0 - (1E0 - magnetic_flux_density_boundary / magnetic_flux_density_end) \
                    * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_cold_temperature_end) \
                        * (magnetic_flux_density_end / magnetic_flux_density_boundary - 1E0))

#    current_density_Knight = current_density_Knight \
#        - electron_hot_number_density_boundary * np.sqrt(electron_hot_temperature_boundary / 2E0 / np.pi / electron_mass) \
#            * (1E0 - (1E0 - magnetic_flux_density_boundary / magnetic_flux_density_end) \
#                * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_hot_temperature_boundary) \
#                    * (magnetic_flux_density_end / magnetic_flux_density_boundary - 1E0))
#
#    current_density_Knight = current_density_Knight \
#        - electron_hot_number_density_end * np.sqrt(electron_hot_temperature_end / 2E0 / np.pi / electron_mass) \
#            * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_hot_temperature_end) \
#                * (1E0 - (1E0 - magnetic_flux_density_boundary / magnetic_flux_density_end) \
#                    * np.exp(- elementary_charge * (electrostatic_potential_end - electrostatic_potential_boundary) / electron_hot_temperature_end) \
#                        * (magnetic_flux_density_end / magnetic_flux_density_boundary - 1E0))

    current_density_Knight = - elementary_charge * magnetic_flux_density_half / magnetic_flux_density_boundary * current_density_Knight

    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'magnetic equator ← MLAT [degree] → North', ylabel=r'Current density (abs) [$\mathrm{\mu A / m^{2}}$]', yscale=r'log') #
    ax.plot(mlat_degree_half, np.abs(current_density_half*1E6), c=r'blue', linestyle=r'solid', linewidth=r'4', label=r'PDS')
    ax.plot(mlat_degree_half, np.abs(current_density_Knight*1E6), c=r'orange', linestyle=r'solid', linewidth=r'4', label=r'Knight relation')
    ax.minorticks_on()
    ax.grid(which=r'both', alpha=0.5)
    ax.legend()

if (channel == 11):
    length = int((1+len(mlat_degree))/2)

    mlat_degree_half = mlat_degree[length-1:]
    length2planet_half = length2planet[length-1:] / 7.1492E7 - 1E0
    electrostatic_potential_half = electrostatic_potential[length-1:]
    number_density_half = number_density[:, length-1:]*1E-6

    electron_perpendicular_pressure = pressure_perp[1, length-1:] + pressure_perp[3, length-1:] + pressure_perp[8, length-1:] + pressure_perp[9, length-1:]
    ion_perpendicular_pressure = pressure_perp[0, length-1:] + pressure_perp[2, length-1:] + pressure_perp[4, length-1:] + pressure_perp[5, length-1:] + pressure_perp[6, length-1:] + pressure_perp[7, length-1:]
    all_perpendicular_pressure = electron_perpendicular_pressure + ion_perpendicular_pressure

    electron_parallel_pressure = pressure_para[1, length-1:] + pressure_para[3, length-1:] + pressure_para[8, length-1:] + pressure_para[9, length-1:]
    ion_parallel_pressure = pressure_para[0, length-1:] + pressure_para[2, length-1:] + pressure_para[4, length-1:] + pressure_para[5, length-1:] + pressure_para[6, length-1:] + pressure_para[7, length-1:]
    all_parallel_pressure = electron_parallel_pressure + ion_parallel_pressure

    electron_total_pressure = electron_parallel_pressure / 3E0 + electron_perpendicular_pressure * 2E0 / 3E0
    ion_total_pressure = ion_parallel_pressure / 3E0 + ion_perpendicular_pressure * 2E0 / 3E0
    all_total_pressure = electron_total_pressure + ion_total_pressure

    electron_para_per_perp_pressure = electron_parallel_pressure / electron_perpendicular_pressure
    ion_para_per_perp_pressure = ion_parallel_pressure / ion_perpendicular_pressure
    all_para_per_perp_pressure = all_parallel_pressure / all_perpendicular_pressure

    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.usetex'] = True
    plt.rcParams["font.size"] = 25

    fig = plt.figure(tight_layout=True)
    gs = fig.add_gridspec(9, 7)
    ax1 = fig.add_subplot(gs[0:3, 0:5], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'number density $n$ [$\mathrm{cm^{-3}}$]', yscale='log', ylim=(1E-2, 1E5))
    ax2 = fig.add_subplot(gs[3:5, 0:], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'electrostatic potential $\phi$ [$\mathrm{kV}$]')
    ax4 = fig.add_subplot(gs[5, 0:], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'$\phi$ [$\mathrm{V}$]')
    ax5 = fig.add_subplot(gs[6:, 0:2], yscale='log', xlabel=r'MLAT [$\mathrm{deg}$]', ylabel=r'pressure [$\mathrm{nPa}$]')
    ax6 = fig.add_subplot(gs[6:, 2:4], yscale='log', xlabel=r'MLAT [$\mathrm{deg}$]')
    ax7 = fig.add_subplot(gs[6:, 4:6], yscale='log', xlabel=r'MLAT [$\mathrm{deg}$]')

    ax1.set_title(r'(a)', x=-0.1, y=0.95)
    ax1.plot(mlat_degree_half, number_density_half[0, :]+number_density_half[2, :], c='dimgrey', label=r'$\mathrm{H^+}$(Jupiter)', linestyle='solid', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree_half, number_density_half[1, :]+number_density_half[3, :], c='blue', label=r'$\mathrm{e^-}$(Jupiter)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[4, :], c='purple', label=r'$\mathrm{H^+}$(Io)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[5, :], c='orange', label=r'$\mathrm{O^+}$(Io)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[6, :], c='green', label=r'$\mathrm{S^+}$(Io)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[7, :], c='lime', label=r'$\mathrm{S^{2+}}$(Io)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[8, :], c='deepskyblue', label=r'cold $\mathrm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_half, number_density_half[9, :], c='hotpink', label=r'hot $\mathrm{e^-}$(Io)', linestyle='dotted', linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both", alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=30)

    ax2.set_title(r'(b)', x=-0.1, y=0.95)
    ax2.plot(mlat_degree_half, electrostatic_potential_half/1E3, linewidth='4', c='blue')
    ax2.minorticks_on()
    ax2.grid(which="both", alpha=0.3)

    ax4.set_title(r'(c)', x=-0.1, y=0.95)
    ax4.plot(mlat_degree_half, electrostatic_potential_half, linewidth='4', c='blue')
    ax4.set_ylim(-20, 5)
    ax4.minorticks_on()
    ax4.grid(which="both", alpha=0.3)

    ax5.set_title(r'(d) perpendicular $P_{\perp}$')
    ax5.plot(mlat_degree_half, all_perpendicular_pressure*1E9, c='purple', label='total', linewidth='4', alpha=0.7)
    ax5.plot(mlat_degree_half, ion_perpendicular_pressure*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax5.plot(mlat_degree_half, electron_perpendicular_pressure*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax5.minorticks_on()
    ax5.grid(which="both", alpha=0.3)

    ax6.set_title(r'(e) parallel $P_{\parallel}$')
    ax6.plot(mlat_degree_half, all_parallel_pressure*1E9, c='purple', label='total', linewidth='4', alpha=0.7)
    ax6.plot(mlat_degree_half, ion_parallel_pressure*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax6.plot(mlat_degree_half, electron_parallel_pressure*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax6.minorticks_on()
    ax6.grid(which="both", alpha=0.3)

    ax7.set_title(r'(f) total $P$')
    ax7.plot(mlat_degree_half, all_total_pressure*1E9, c='purple', label='total', linewidth='4', alpha=0.7)
    ax7.plot(mlat_degree_half, ion_total_pressure*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax7.plot(mlat_degree_half, electron_total_pressure*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax7.minorticks_on()
    ax7.grid(which="both", alpha=0.3)
    ax7.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=30)

    #ax7.set_title(r'(f) anisotropy $P_{\parallel} / P_{\perp}$')
    #ax7.plot(mlat_degree_half, all_para_per_perp_pressure, c='purple', label='total', linewidth='4', alpha=0.7)
    #ax7.plot(mlat_degree_half, ion_para_per_perp_pressure, c='orange', label='ion', linewidth='4', alpha=0.7)
    #ax7.plot(mlat_degree_half, electron_para_per_perp_pressure, c='blue', label='electron', linewidth='4', alpha=0.7)
    #ax7.minorticks_on()
    #ax7.grid(which="both", alpha=0.3)
    #ax7.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=30)

    plt.subplots_adjust(wspace=1, hspace=1)

if (channel == 12):
    length = int((1+len(mlat_degree))/2)

    coordinate_FA_half = coordinate_FA[length-1:]
    Alfven_speed_half = Alfven_speed[length-1:]

    transit_time = np.zeros(length)
    length2planet_half = length2planet[length-1:] / 7.1492E7 - 1E0
    mlat_degree_half = mlat_degree[length-1:]

    for length_i in range(length-1):
        distance_coordinate_FA = coordinate_FA_half[length_i+1] - coordinate_FA_half[length_i]
        average_Alfven_speed = (Alfven_speed_half[length_i+1] + Alfven_speed_half[length_i]) / 2E0
        transit_time[length_i+1] = transit_time[length_i] + distance_coordinate_FA / average_Alfven_speed

    one_length = np.ones(length)

    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=r'propagation time [$\mathrm{s}$]' ,ylabel=r'MLAT [$\mathrm{deg}$]')
    ax.plot(transit_time, mlat_degree_half, c='blue', linestyle='solid', linewidth='4', label='propagation')
    ax.plot(one_length * max(transit_time), mlat_degree_half, c='black', linestyle='dotted', linewidth='4', alpha=0.5, label=r'time = ' + str(round(max(transit_time), 2)) + r' $s$')
    ax.plot(transit_time, one_length * max(mlat_degree_half), c='black', linestyle='dotted', linewidth='4', alpha=0.5, label=r'MLAT = ' + str(round(max(mlat_degree_half), 2)) + r' degree')
    ax.minorticks_on()
    ax.grid(which="both")
    ax.legend()

if (channel == 13):
    length = int((1+len(mlat_degree))/2)

    mlat_degree_half = mlat_degree[length-1:]
    length2planet_half = length2planet[length-1:] / 7.1492E7 - 1E0

    electron_perpendicular_pressure = pressure_perp[1, length-1:] + pressure_perp[3, length-1:] + pressure_perp[8, length-1:] + pressure_perp[9, length-1:]
    ion_perpendicular_pressure = pressure_perp[0, length-1:] + pressure_perp[2, length-1:] + pressure_perp[4, length-1:] + pressure_perp[5, length-1:] + pressure_perp[6, length-1:] + pressure_perp[7, length-1:]
    all_perpendicular_pressure = electron_perpendicular_pressure + ion_perpendicular_pressure

    electron_parallel_pressure = pressure_para[1, length-1:] + pressure_para[3, length-1:] + pressure_para[8, length-1:] + pressure_para[9, length-1:]
    ion_parallel_pressure = pressure_para[0, length-1:] + pressure_para[2, length-1:] + pressure_para[4, length-1:] + pressure_para[5, length-1:] + pressure_para[6, length-1:] + pressure_para[7, length-1:]
    all_parallel_pressure = electron_parallel_pressure + ion_parallel_pressure

    electron_total_pressure = electron_parallel_pressure / 3E0 + electron_perpendicular_pressure * 2E0 / 3E0
    ion_total_pressure = ion_parallel_pressure / 3E0 + ion_perpendicular_pressure * 2E0 / 3E0
    all_total_pressure = electron_total_pressure + ion_total_pressure

    magnetic_constant = 1.25663706212E-6
    magnetic_pressure = magnetic_flux_density[length-1:]**2E0 / 2E0 / magnetic_constant

    number_density_half = number_density[:, length-1:]
    ion_mass_density = (number_density_half[0, :] + number_density_half[2, :] + number_density_half[4, :]) * 1.67262192369E-27 + number_density_half[5, :] * 2.677950266103E-26 + number_density_half[6, :] * 5.355991626103E-26 + number_density_half[7, :] * 5.355900532206E-26
    electron_mass_density = (number_density_half[1, :] + number_density_half[3, :] + number_density_half[8, :] + number_density_half[9, :]) * 9.1093837015E-31
    electron2ion_mass_ratio =electron_mass_density / ion_mass_density

    plt.rcParams["font.size"] = 60
    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.set_title(r'perpendicular')
    ax1.set_xlabel(r'MLAT [$\mathrm{deg}$]')
    ax1.set_ylabel(r'plasma β')
    ax1.set_yscale(r'log')
    ax1.plot(mlat_degree_half, electron2ion_mass_ratio, c='dimgrey', label=r'$m_{e}/m_{i}$', linestyle='dotted', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree_half, electron_perpendicular_pressure / magnetic_pressure, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree_half, ion_perpendicular_pressure / magnetic_pressure, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax1.plot(mlat_degree_half, all_perpendicular_pressure / magnetic_pressure, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax1.minorticks_on()
    ax1.grid(which="both")

    ax2.set_title(r'parallel')
    ax2.set_xlabel(r'MLAT [$\mathrm{deg}$]')
    ax2.set_yscale(r'log')
    ax2.plot(mlat_degree_half, electron2ion_mass_ratio, c='dimgrey', label=r'$m_{e}/m_{i}$', linestyle='dotted', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree_half, electron_parallel_pressure / magnetic_pressure, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree_half, ion_parallel_pressure / magnetic_pressure, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax2.plot(mlat_degree_half, all_parallel_pressure / magnetic_pressure, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax2.minorticks_on()
    ax2.grid(which="both")

    ax3.set_title(r'total')
    ax3.set_xlabel(r'MLAT [$\mathrm{deg}$]')
    ax3.set_yscale(r'log')
    ax3.plot(mlat_degree_half, electron2ion_mass_ratio, c='dimgrey', label=r'$m_{e}/m_{i}$', linestyle='dotted', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree_half, electron_total_pressure / magnetic_pressure, c='blue', label=r'electron', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree_half, ion_total_pressure / magnetic_pressure, c='orange', label=r'ion', linewidth='4', alpha=0.7)
    ax3.plot(mlat_degree_half, all_total_pressure / magnetic_pressure, c='purple', label=r'all', linewidth='4', alpha=0.7)
    ax3.minorticks_on()
    ax3.grid(which="both")
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    

    plt.subplots_adjust(wspace=0.4, hspace=0.6)


plt.tight_layout()
plt.show()