import numpy as np
import matplotlib.pyplot as plt

data_Case_1 = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_050_BC_001_min_088.csv", delimiter=',', unpack=True)
data_Case_2 = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_050_BC_004_min_074.csv", delimiter=',', unpack=True)
data_Case_3 = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_027_BC_001_min_073.csv", delimiter=',', unpack=True)
data_Case_4 = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_027_BC_004_min_063.csv", delimiter=',', unpack=True)

series_number = 10

magnetic_constant = 1.25663706212E-6

#data_Case_1
dummy_1 = data_Case_1[3, :]
half_length_Case_1 = int((1+len(dummy_1))/2)
mlat_degree_Case_1 = data_Case_1[3, half_length_Case_1-1:]

perpendicular_pressure_Case_1 = data_Case_1[3 * series_number + 10 : 4 * series_number + 10, half_length_Case_1-1:]
electron_perpendicular_pressure_Case_1 = perpendicular_pressure_Case_1[1, :] + perpendicular_pressure_Case_1[3, :] + perpendicular_pressure_Case_1[8, :] + perpendicular_pressure_Case_1[9, :]
ion_perpendicular_pressure_Case_1 = perpendicular_pressure_Case_1[0, :] + perpendicular_pressure_Case_1[2, :] + perpendicular_pressure_Case_1[4, :] + perpendicular_pressure_Case_1[5, :] + perpendicular_pressure_Case_1[6, :] + perpendicular_pressure_Case_1[7, :]
all_perpendicular_pressure_Case_1 = electron_perpendicular_pressure_Case_1 + ion_perpendicular_pressure_Case_1

parallel_pressure_Case_1 = data_Case_1[4 * series_number + 10 : 5 * series_number + 10, half_length_Case_1-1:]
electron_parallel_pressure_Case_1 = parallel_pressure_Case_1[1, :] + parallel_pressure_Case_1[3, :] + parallel_pressure_Case_1[8, :] + parallel_pressure_Case_1[9, :]
ion_parallel_pressure_Case_1 = parallel_pressure_Case_1[0, :] + parallel_pressure_Case_1[2, :] + parallel_pressure_Case_1[4, :] + parallel_pressure_Case_1[5, :] + parallel_pressure_Case_1[6, :] + parallel_pressure_Case_1[7, :]
all_parallel_pressure_Case_1 = electron_parallel_pressure_Case_1 + ion_parallel_pressure_Case_1

electron_total_pressure_Case_1 = electron_parallel_pressure_Case_1 / 3E0 + electron_perpendicular_pressure_Case_1 * 2E0 / 3E0
ion_total_pressure_Case_1 = ion_parallel_pressure_Case_1 / 3E0 + ion_perpendicular_pressure_Case_1 * 2E0 / 3E0
all_total_pressure_Case_1 = all_parallel_pressure_Case_1 / 3E0 + all_perpendicular_pressure_Case_1 * 2E0 / 3E0

magnetic_flux_density_Case_1 = data_Case_1[4, half_length_Case_1-1:]
magnetic_pressure_Case_1 = magnetic_flux_density_Case_1**2E0 / 2E0 / magnetic_constant

electron_total_beta_Case_1 = electron_total_pressure_Case_1 / magnetic_pressure_Case_1
ion_total_beta_Case_1 = ion_total_pressure_Case_1 / magnetic_pressure_Case_1
all_total_beta_Case_1 = all_total_pressure_Case_1 / magnetic_pressure_Case_1

number_density_Case_1 = data_Case_1[7:series_number + 7, half_length_Case_1-1:]
ion_mass_density_Case_1 = (number_density_Case_1[0, :] + number_density_Case_1[2, :] + number_density_Case_1[4, :]) * 1.67262192369E-27 + number_density_Case_1[5, :] * 2.677950266103E-26 + number_density_Case_1[6, :] * 5.355991626103E-26 + number_density_Case_1[7, :] * 5.355900532206E-26
electron_mass_density_Case_1 = (number_density_Case_1[1, :] + number_density_Case_1[3, :] + number_density_Case_1[8, :] + number_density_Case_1[9, :]) * 9.1093837015E-31
electron2ion_mass_ratio_Case_1 =electron_mass_density_Case_1 / ion_mass_density_Case_1


#data_Case_2
dummy_2 = data_Case_2[3, :]
half_length_Case_2 = int((1+len(dummy_2))/2)
mlat_degree_Case_2 = data_Case_2[3, half_length_Case_2-1:]

perpendicular_pressure_Case_2 = data_Case_2[3 * series_number + 10 : 4 * series_number + 10, half_length_Case_2-1:]
electron_perpendicular_pressure_Case_2 = perpendicular_pressure_Case_2[1, :] + perpendicular_pressure_Case_2[3, :] + perpendicular_pressure_Case_2[8, :] + perpendicular_pressure_Case_2[9, :]
ion_perpendicular_pressure_Case_2 = perpendicular_pressure_Case_2[0, :] + perpendicular_pressure_Case_2[2, :] + perpendicular_pressure_Case_2[4, :] + perpendicular_pressure_Case_2[5, :] + perpendicular_pressure_Case_2[6, :] + perpendicular_pressure_Case_2[7, :]
all_perpendicular_pressure_Case_2 = electron_perpendicular_pressure_Case_2 + ion_perpendicular_pressure_Case_2

parallel_pressure_Case_2 = data_Case_2[4 * series_number + 10 : 5 * series_number + 10, half_length_Case_2-1:]
electron_parallel_pressure_Case_2 = parallel_pressure_Case_2[1, :] + parallel_pressure_Case_2[3, :] + parallel_pressure_Case_2[8, :] + parallel_pressure_Case_2[9, :]
ion_parallel_pressure_Case_2 = parallel_pressure_Case_2[0, :] + parallel_pressure_Case_2[2, :] + parallel_pressure_Case_2[4, :] + parallel_pressure_Case_2[5, :] + parallel_pressure_Case_2[6, :] + parallel_pressure_Case_2[7, :]
all_parallel_pressure_Case_2 = electron_parallel_pressure_Case_2 + ion_parallel_pressure_Case_2

electron_total_pressure_Case_2 = electron_parallel_pressure_Case_2 / 3E0 + electron_perpendicular_pressure_Case_2 * 2E0 / 3E0
ion_total_pressure_Case_2 = ion_parallel_pressure_Case_2 / 3E0 + ion_perpendicular_pressure_Case_2 * 2E0 / 3E0
all_total_pressure_Case_2 = all_parallel_pressure_Case_2 / 3E0 + all_perpendicular_pressure_Case_2 * 2E0 / 3E0

magnetic_flux_density_Case_2 = data_Case_2[4, half_length_Case_2-1:]
magnetic_pressure_Case_2 = magnetic_flux_density_Case_2**2E0 / 2E0 / magnetic_constant

electron_total_beta_Case_2 = electron_total_pressure_Case_2 / magnetic_pressure_Case_2
ion_total_beta_Case_2 = ion_total_pressure_Case_2 / magnetic_pressure_Case_2
all_total_beta_Case_2 = all_total_pressure_Case_2 / magnetic_pressure_Case_2

number_density_Case_2 = data_Case_2[7:series_number + 7, half_length_Case_2-1:]
ion_mass_density_Case_2 = (number_density_Case_2[0, :] + number_density_Case_2[2, :] + number_density_Case_2[4, :]) * 1.67262192369E-27 + number_density_Case_2[5, :] * 2.677950266103E-26 + number_density_Case_2[6, :] * 5.355991626103E-26 + number_density_Case_2[7, :] * 5.355900532206E-26
electron_mass_density_Case_2 = (number_density_Case_2[1, :] + number_density_Case_2[3, :] + number_density_Case_2[8, :] + number_density_Case_2[9, :]) * 9.1093837015E-31
electron2ion_mass_ratio_Case_2 =electron_mass_density_Case_2 / ion_mass_density_Case_2


#data_Case_3
dummy_3 = data_Case_3[3, :]
half_length_Case_3 = int((1+len(dummy_3))/2)
mlat_degree_Case_3 = data_Case_3[3, half_length_Case_3-1:]

perpendicular_pressure_Case_3 = data_Case_3[3 * series_number + 10 : 4 * series_number + 10, half_length_Case_3-1:]
electron_perpendicular_pressure_Case_3 = perpendicular_pressure_Case_3[1, :] + perpendicular_pressure_Case_3[3, :] + perpendicular_pressure_Case_3[8, :] + perpendicular_pressure_Case_3[9, :]
ion_perpendicular_pressure_Case_3 = perpendicular_pressure_Case_3[0, :] + perpendicular_pressure_Case_3[2, :] + perpendicular_pressure_Case_3[4, :] + perpendicular_pressure_Case_3[5, :] + perpendicular_pressure_Case_3[6, :] + perpendicular_pressure_Case_3[7, :]
all_perpendicular_pressure_Case_3 = electron_perpendicular_pressure_Case_3 + ion_perpendicular_pressure_Case_3

parallel_pressure_Case_3 = data_Case_3[4 * series_number + 10 : 5 * series_number + 10, half_length_Case_3-1:]
electron_parallel_pressure_Case_3 = parallel_pressure_Case_3[1, :] + parallel_pressure_Case_3[3, :] + parallel_pressure_Case_3[8, :] + parallel_pressure_Case_3[9, :]
ion_parallel_pressure_Case_3 = parallel_pressure_Case_3[0, :] + parallel_pressure_Case_3[2, :] + parallel_pressure_Case_3[4, :] + parallel_pressure_Case_3[5, :] + parallel_pressure_Case_3[6, :] + parallel_pressure_Case_3[7, :]
all_parallel_pressure_Case_3 = electron_parallel_pressure_Case_3 + ion_parallel_pressure_Case_3

electron_total_pressure_Case_3 = electron_parallel_pressure_Case_3 / 3E0 + electron_perpendicular_pressure_Case_3 * 2E0 / 3E0
ion_total_pressure_Case_3 = ion_parallel_pressure_Case_3 / 3E0 + ion_perpendicular_pressure_Case_3 * 2E0 / 3E0
all_total_pressure_Case_3 = all_parallel_pressure_Case_3 / 3E0 + all_perpendicular_pressure_Case_3 * 2E0 / 3E0

magnetic_flux_density_Case_3 = data_Case_3[4, half_length_Case_3-1:]
magnetic_pressure_Case_3 = magnetic_flux_density_Case_3**2E0 / 2E0 / magnetic_constant

electron_total_beta_Case_3 = electron_total_pressure_Case_3 / magnetic_pressure_Case_3
ion_total_beta_Case_3 = ion_total_pressure_Case_3 / magnetic_pressure_Case_3
all_total_beta_Case_3 = all_total_pressure_Case_3 / magnetic_pressure_Case_3

number_density_Case_3 = data_Case_3[7:series_number + 7, half_length_Case_3-1:]
ion_mass_density_Case_3 = (number_density_Case_3[0, :] + number_density_Case_3[2, :] + number_density_Case_3[4, :]) * 1.67262192369E-27 + number_density_Case_3[5, :] * 2.677950266103E-26 + number_density_Case_3[6, :] * 5.355991626103E-26 + number_density_Case_3[7, :] * 5.355900532206E-26
electron_mass_density_Case_3 = (number_density_Case_3[1, :] + number_density_Case_3[3, :] + number_density_Case_3[8, :] + number_density_Case_3[9, :]) * 9.1093837015E-31
electron2ion_mass_ratio_Case_3 =electron_mass_density_Case_3 / ion_mass_density_Case_3


#data_Case_4
dummy_4 = data_Case_4[3, :]
half_length_Case_4 = int((1+len(dummy_4))/2)
mlat_degree_Case_4 = data_Case_4[3, half_length_Case_4-1:]

perpendicular_pressure_Case_4 = data_Case_4[3 * series_number + 10 : 4 * series_number + 10, half_length_Case_4-1:]
electron_perpendicular_pressure_Case_4 = perpendicular_pressure_Case_4[1, :] + perpendicular_pressure_Case_4[3, :] + perpendicular_pressure_Case_4[8, :] + perpendicular_pressure_Case_4[9, :]
ion_perpendicular_pressure_Case_4 = perpendicular_pressure_Case_4[0, :] + perpendicular_pressure_Case_4[2, :] + perpendicular_pressure_Case_4[4, :] + perpendicular_pressure_Case_4[5, :] + perpendicular_pressure_Case_4[6, :] + perpendicular_pressure_Case_4[7, :]
all_perpendicular_pressure_Case_4 = electron_perpendicular_pressure_Case_4 + ion_perpendicular_pressure_Case_4

parallel_pressure_Case_4 = data_Case_4[4 * series_number + 10 : 5 * series_number + 10, half_length_Case_4-1:]
electron_parallel_pressure_Case_4 = parallel_pressure_Case_4[1, :] + parallel_pressure_Case_4[3, :] + parallel_pressure_Case_4[8, :] + parallel_pressure_Case_4[9, :]
ion_parallel_pressure_Case_4 = parallel_pressure_Case_4[0, :] + parallel_pressure_Case_4[2, :] + parallel_pressure_Case_4[4, :] + parallel_pressure_Case_4[5, :] + parallel_pressure_Case_4[6, :] + parallel_pressure_Case_4[7, :]
all_parallel_pressure_Case_4 = electron_parallel_pressure_Case_4 + ion_parallel_pressure_Case_4

electron_total_pressure_Case_4 = electron_parallel_pressure_Case_4 / 3E0 + electron_perpendicular_pressure_Case_4 * 2E0 / 3E0
ion_total_pressure_Case_4 = ion_parallel_pressure_Case_4 / 3E0 + ion_perpendicular_pressure_Case_4 * 2E0 / 3E0
all_total_pressure_Case_4 = all_parallel_pressure_Case_4 / 3E0 + all_perpendicular_pressure_Case_4 * 2E0 / 3E0

magnetic_flux_density_Case_4 = data_Case_4[4, half_length_Case_4-1:]
magnetic_pressure_Case_4 = magnetic_flux_density_Case_4**2E0 / 2E0 / magnetic_constant

electron_total_beta_Case_4 = electron_total_pressure_Case_4 / magnetic_pressure_Case_4
ion_total_beta_Case_4 = ion_total_pressure_Case_4 / magnetic_pressure_Case_4
all_total_beta_Case_4 = all_total_pressure_Case_4 / magnetic_pressure_Case_4

number_density_Case_4 = data_Case_4[7:series_number + 7, half_length_Case_4-1:]
ion_mass_density_Case_4 = (number_density_Case_4[0, :] + number_density_Case_4[2, :] + number_density_Case_4[4, :]) * 1.67262192369E-27 + number_density_Case_4[5, :] * 2.677950266103E-26 + number_density_Case_4[6, :] * 5.355991626103E-26 + number_density_Case_4[7, :] * 5.355900532206E-26
electron_mass_density_Case_4 = (number_density_Case_4[1, :] + number_density_Case_4[3, :] + number_density_Case_4[8, :] + number_density_Case_4[9, :]) * 9.1093837015E-31
electron2ion_mass_ratio_Case_4 =electron_mass_density_Case_4 / ion_mass_density_Case_4


#plot
plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 25

fig = plt.figure(tight_layout=True)
gs = fig.add_gridspec(2, 2)
ax1 = fig.add_subplot(gs[0, 0], xlabel=r'MLAT [$\mathrm{deg}$]', ylabel=r'plasma $\beta$', yscale='log', title=r'(a) Case 1')
ax2 = fig.add_subplot(gs[0, 1], xlabel=r'MLAT [$\mathrm{deg}$]', ylabel=r'plasma $\beta$', yscale='log', title=r'(b) Case 2')
ax3 = fig.add_subplot(gs[1, 0], xlabel=r'MLAT [$\mathrm{deg}$]', ylabel=r'plasma $\beta$', yscale='log', title=r'(c) Case 3')
ax4 = fig.add_subplot(gs[1, 1], xlabel=r'MLAT [$\mathrm{deg}$]', ylabel=r'plasma $\beta$', yscale='log', title=r'(d) Case 4')

ax1.plot(mlat_degree_Case_1, electron2ion_mass_ratio_Case_1, c=r'dimgrey', linestyle=r'dotted', linewidth=r'4', alpha=0.7, label=r'$m_{e}/m_{i}$')
ax1.plot(mlat_degree_Case_1, electron_total_beta_Case_1, c=r'blue', linewidth=r'4', alpha=0.7, label=r'electron')
ax1.plot(mlat_degree_Case_1, ion_total_beta_Case_1, c=r'orange', linewidth=r'4', alpha=0.7, label=r'ion')
ax1.plot(mlat_degree_Case_1, all_total_beta_Case_1, c=r'purple', linewidth=r'4', alpha=0.7, label=r'all')
ax1.minorticks_on()
ax1.grid(which="both", alpha=0.5)
ax1.set_yticks(np.logspace(-11, -1, 11))
ax1.legend()

ax2.plot(mlat_degree_Case_2, electron2ion_mass_ratio_Case_2, c=r'dimgrey', linestyle=r'dotted', linewidth=r'4', alpha=0.7, label=r'$m_{e}/m_{i}$')
ax2.plot(mlat_degree_Case_2, electron_total_beta_Case_2, c=r'blue', linewidth=r'4', alpha=0.7, label=r'electron')
ax2.plot(mlat_degree_Case_2, ion_total_beta_Case_2, c=r'orange', linewidth=r'4', alpha=0.7, label=r'ion')
ax2.plot(mlat_degree_Case_2, all_total_beta_Case_2, c=r'purple', linewidth=r'4', alpha=0.7, label=r'all')
ax2.minorticks_on()
ax2.grid(which="both", alpha=0.5)
ax2.set_yticks(np.logspace(-11, -1, 11))
ax2.legend()

ax3.plot(mlat_degree_Case_3, electron2ion_mass_ratio_Case_3, c=r'dimgrey', linestyle=r'dotted', linewidth=r'4', alpha=0.7, label=r'$m_{e}/m_{i}$')
ax3.plot(mlat_degree_Case_3, electron_total_beta_Case_3, c=r'blue', linewidth=r'4', alpha=0.7, label=r'electron')
ax3.plot(mlat_degree_Case_3, ion_total_beta_Case_3, c=r'orange', linewidth=r'4', alpha=0.7, label=r'ion')
ax3.plot(mlat_degree_Case_3, all_total_beta_Case_3, c=r'purple', linewidth=r'4', alpha=0.7, label=r'all')
ax3.minorticks_on()
ax3.grid(which="both", alpha=0.5)
ax3.set_yticks(np.logspace(-11, -1, 11))
ax3.legend()

ax4.plot(mlat_degree_Case_4, electron2ion_mass_ratio_Case_4, c=r'dimgrey', linestyle=r'dotted', linewidth=r'4', alpha=0.7, label=r'$m_{e}/m_{i}$')
ax4.plot(mlat_degree_Case_4, electron_total_beta_Case_4, c=r'blue', linewidth=r'4', alpha=0.7, label=r'electron')
ax4.plot(mlat_degree_Case_4, ion_total_beta_Case_4, c=r'orange', linewidth=r'4', alpha=0.7, label=r'ion')
ax4.plot(mlat_degree_Case_4, all_total_beta_Case_4, c=r'purple', linewidth=r'4', alpha=0.7, label=r'all')
ax4.minorticks_on()
ax4.grid(which="both", alpha=0.5)
ax4.set_yticks(np.logspace(-11, -1, 11))
ax4.legend()

plt.tight_layout()
plt.show()