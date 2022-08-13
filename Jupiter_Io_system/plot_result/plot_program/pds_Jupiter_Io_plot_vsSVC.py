import numpy as np
import matplotlib.pyplot as plt

data_pds = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_050_BC_001_min_088.csv", delimiter=',', unpack=True)

data_pds_2 = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results_all/result_all_027_050_BC_004_min_074.csv", delimiter=',', unpack=True)

data_svc = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/results/result_number_density_027_050_BC_001_min_070_SVC_represent.csv", delimiter=',', unpack=True, skip_header=1, skip_footer=10)

data_svc_Matsuda_numberdensity = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/SVC_J/SVC_penum_39_50_75.csv", delimiter=',', unpack=True)
data_svc_Matsuda_initialcondition = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/SVC_J/SVC_IC_J_30kV_75_50.csv", delimiter=',', unpack=True)
data_svc_Matsuda_magneticfluxdensity = np.genfromtxt(r"/home/satanka/Documents/fort/pds_kai/pds_J_kai/SVC_J/mag_FA.csv", delimiter=',', unpack=True)


channel = 2
#1: for paper (number density, electrostatic potential), 2: for paper (svc_Matsuda)

series_number = 10

mlat_degree_pds = data_pds[3, :]
electrostatic_potential_pds = data_pds[6, :]
numberdensity_pds = data_pds[7:series_number + 7, :]/1E6
numberdensity_H_pds = numberdensity_pds[0, :] + numberdensity_pds[2, :] + numberdensity_pds[4, :]
numberdensity_O_pds = numberdensity_pds[5, :]
numberdensity_e_pds = numberdensity_pds[1, :] + numberdensity_pds[3, :] + numberdensity_pds[8, :] + numberdensity_pds[9, :]

mlat_degree_pds_2 = data_pds_2[3, :]
electrostatic_potential_pds_2 = data_pds_2[6, :]
numberdensity_pds_2 = data_pds_2[7:series_number + 7, :]/1E6
numberdensity_H_pds_2 = numberdensity_pds_2[0, :] + numberdensity_pds_2[2, :] + numberdensity_pds_2[4, :]
numberdensity_O_pds_2 = numberdensity_pds_2[5, :]
numberdensity_e_pds_2 = numberdensity_pds_2[1, :] + numberdensity_pds_2[3, :] + numberdensity_pds_2[8, :] + numberdensity_pds_2[9, :]

mlat_degree_svc = data_svc[3, :]
electrostatic_potential_svc = data_svc[6, :]
numberdensity_svc = data_svc[7:series_number + 7, :]/1E6
numberdensity_H_svc = numberdensity_svc[0, :] + numberdensity_svc[2, :] + numberdensity_svc[4, :]
numberdensity_O_svc = numberdensity_svc[5, :]
numberdensity_e_svc = numberdensity_svc[1, :] + numberdensity_svc[3, :] + numberdensity_svc[8, :] + numberdensity_svc[9, :]

coordinate_FA_svc_Matsuda = data_svc_Matsuda_initialcondition[2, 1:]
planet_radius = 7.1492E7
planet_l_shell = 5.8476
radius_equator = planet_radius * planet_l_shell
length_svc_Matsuda = len(coordinate_FA_svc_Matsuda)
mlat_svc_Matsuda = np.zeros(length_svc_Matsuda)
electrostatic_potential_svc_Matsuda = data_svc_Matsuda_numberdensity[0, :]
numberdensity_svc_Matsuda = data_svc_Matsuda_numberdensity[1:9, :]/1E6
numberdensity_H_svc_Matsuda = numberdensity_svc_Matsuda[0, :] + numberdensity_svc_Matsuda[5, :]
numberdensity_O_svc_Matsuda = numberdensity_svc_Matsuda[2, :]
numberdensity_e_svc_Matsuda = numberdensity_svc_Matsuda[1, :] + numberdensity_svc_Matsuda[6, :] + numberdensity_svc_Matsuda[7, :]

#MLAT calculation
for grid_i in range(length_svc_Matsuda):
    mlat_0 = 1.
    for iteration_i in range(1000000):
        if (iteration_i == 1000000):
            print("Error!: solution is not found. z_position = " + str(coordinate_FA_svc_Matsuda[grid_i]))

        ff = radius_equator * (5E-1 * np.sin(mlat_0) * np.sqrt(3E0 * np.sin(mlat_0)**2E0 + 1E0) + 1E0 / 2E0 / np.sqrt(3E0) * np.arcsinh(np.sqrt(3E0) * np.sin(mlat_0))) - coordinate_FA_svc_Matsuda[grid_i]
        gg = radius_equator * np.cos(mlat_0) * np.sqrt(3E0 * np.sin(mlat_0)**2E0 + 1E0)

        mlat_1 = float(mlat_0 - ff / gg)

        if (abs(mlat_1 - mlat_0) < 1E-5):
            break

        mlat_0 = mlat_1

    mlat_svc_Matsuda[grid_i] = mlat_1

mlat_degree_svc_Matsuda = mlat_svc_Matsuda / np.pi * 180E0


if (channel == 1):
    half_length_pds = int((1+len(mlat_degree_pds))/2)
    half_length_pds_2 = int((1+len(mlat_degree_pds_2))/2)
    half_length_svc = int((1+len(mlat_degree_svc))/2)

    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.usetex'] = True
    plt.rcParams["font.size"] = 25

    fig = plt.figure(tight_layout=True)
    gs = fig.add_gridspec(4, 1)
    ax1 = fig.add_subplot(gs[0:2, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'number density $n$ [$\rm{cm^{-3}}$]', yscale='log', ylim=(1.E-2, 1.E5))
    ax2 = fig.add_subplot(gs[2, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'electrostatic potential $\phi$ [$\mathrm{kV}$]')
    ax3 = fig.add_subplot(gs[3, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'$\phi$ [$\mathrm{V}$]')

    ax1.set_title(r'(a)', x=-0.1, y=0.95)
    #ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_H_pds[half_length_pds-1:], c='red', label=r'$\rm{H^+}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_e_pds[half_length_pds-1:], c='blue', label=r'$\rm{e^-}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_O_pds[half_length_pds-1:], c='orange', label=r'$\rm{O^+}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_H_pds_2[half_length_pds_2-1:], c='red', label=r'$\rm{H^+}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_e_pds_2[half_length_pds_2-1:], c='blue', label=r'$\rm{e^-}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_O_pds_2[half_length_pds_2-1:], c='orange', label=r'$\rm{O^+}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    ax1.plot(mlat_degree_svc[half_length_svc-1:], numberdensity_H_svc[half_length_svc-1:], c='red', label=r'$\rm{H^+}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_svc[half_length_svc-1:], numberdensity_e_svc[half_length_svc-1:], c='blue', label=r'$\rm{e^-}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_svc[half_length_svc-1:], numberdensity_O_svc[half_length_svc-1:], c='orange', label=r'$\rm{O^+}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both", alpha=0.3)
    ax1.legend()

    ax2.set_title(r'(b)', x=-0.1, y=0.95)
    #ax2.plot(mlat_degree_pds[half_length_pds-1:], electrostatic_potential_pds[half_length_pds-1:]/1E3, linewidth='4', linestyle='solid', c='red', label='PDS Case 1', alpha=0.5)
    #ax2.plot(mlat_degree_pds_2[half_length_pds_2-1:], electrostatic_potential_pds_2[half_length_pds_2-1:]/1E3, linewidth='4', linestyle='-.', c='orange', label='PDS Case 2', alpha=0.5)
    ax2.plot(mlat_degree_svc[half_length_svc-1:], electrostatic_potential_svc[half_length_svc-1:]/1E3, linewidth='4', linestyle='dotted', c='blue', label='SVC')
    ax2.minorticks_on()
    ax2.grid(which="both", alpha=0.3)
    ax2.legend()

    ax3.set_title(r'(c)', x=-0.1, y=0.95)
    ax3.set_ylim(-20, 5)
    #ax3.plot(mlat_degree_pds[half_length_pds-1:], electrostatic_potential_pds[half_length_pds-1:], linewidth='4', linestyle='solid', c='red', label='PDS Case 1', alpha=0.5)
    #ax3.plot(mlat_degree_pds_2[half_length_pds_2-1:], electrostatic_potential_pds_2[half_length_pds_2-1:], linewidth='4', linestyle='-.', c='orange', label='PDS Case 2', alpha=0.5)
    ax3.plot(mlat_degree_svc[half_length_svc-1:], electrostatic_potential_svc[half_length_svc-1:], linewidth='4', linestyle='dotted', c='blue', label='SVC')
    ax3.minorticks_on()
    ax3.grid(which="both", alpha=0.3)
    ax3.legend()

    plt.tight_layout()
    plt.show()


if (channel == 2):
    half_length_pds = int((1+len(mlat_degree_pds))/2)
    half_length_pds_2 = int((1+len(mlat_degree_pds_2))/2)
    half_length_svc = int((1+len(mlat_degree_svc))/2)

    plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
    plt.rcParams['text.usetex'] = True
    plt.rcParams["font.size"] = 25

    fig = plt.figure(tight_layout=True)
    gs = fig.add_gridspec(4, 1)
    ax1 = fig.add_subplot(gs[0:2, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'number density $n$ [$\rm{cm^{-3}}$]', yscale='log', ylim=(1.E-2, 1.E5))
    ax2 = fig.add_subplot(gs[2, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'electrostatic potential $\phi$ [$\mathrm{kV}$]')
    ax3 = fig.add_subplot(gs[3, 0], xlabel=r'magnetic equator ← MLAT [$\mathrm{deg}$] → North', ylabel=r'$\phi$ [$\mathrm{V}$]')

    ax1.set_title(r'(a)', x=-0.1, y=0.95)
    ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_H_pds[half_length_pds-1:], c='red', label=r'$\rm{H^+}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_e_pds[half_length_pds-1:], c='blue', label=r'$\rm{e^-}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    ax1.plot(mlat_degree_pds[half_length_pds-1:], numberdensity_O_pds[half_length_pds-1:], c='orange', label=r'$\rm{O^+}$ (PDS Case 1)', linestyle='solid', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_H_pds_2[half_length_pds_2-1:], c='red', label=r'$\rm{H^+}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_e_pds_2[half_length_pds_2-1:], c='blue', label=r'$\rm{e^-}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    #ax1.plot(mlat_degree_pds_2[half_length_pds_2-1:], numberdensity_O_pds_2[half_length_pds_2-1:], c='orange', label=r'$\rm{O^+}$ (PDS Case 2)', linestyle='-.', linewidth='4', alpha=0.5)
    ax1.plot(mlat_degree_svc_Matsuda, numberdensity_H_svc_Matsuda, c='red', label=r'$\rm{H^+}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_svc_Matsuda, numberdensity_e_svc_Matsuda, c='blue', label=r'$\rm{e^-}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.plot(mlat_degree_svc_Matsuda, numberdensity_O_svc_Matsuda, c='orange', label=r'$\rm{O^+}$ (SVC)', linestyle='dotted', linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which="both", alpha=0.3)
    ax1.legend()

    ax2.set_title(r'(b)', x=-0.1, y=0.95)
    ax2.plot(mlat_degree_pds[half_length_pds-1:], electrostatic_potential_pds[half_length_pds-1:]/1E3, linewidth='4', linestyle='solid', c='red', label='PDS Case 1', alpha=0.5)
    ax2.plot(mlat_degree_pds_2[half_length_pds_2-1:], electrostatic_potential_pds_2[half_length_pds_2-1:]/1E3, linewidth='4', linestyle='-.', c='orange', label='PDS Case 2', alpha=0.5)
    ax2.plot(mlat_degree_svc_Matsuda, electrostatic_potential_svc_Matsuda/1E3, linewidth='4', linestyle='dotted', c='blue', label='SVC')
    ax2.minorticks_on()
    ax2.grid(which="both", alpha=0.3)
    ax2.legend()

    ax3.set_title(r'(c)', x=-0.1, y=0.95)
    ax3.set_ylim(-20, 5)
    ax3.plot(mlat_degree_pds[half_length_pds-1:], electrostatic_potential_pds[half_length_pds-1:], linewidth='4', linestyle='solid', c='red', label='PDS Case 1', alpha=0.5)
    ax3.plot(mlat_degree_pds_2[half_length_pds_2-1:], electrostatic_potential_pds_2[half_length_pds_2-1:], linewidth='4', linestyle='-.', c='orange', label='PDS Case 2', alpha=0.5)
    ax3.plot(mlat_degree_svc_Matsuda, electrostatic_potential_svc_Matsuda, linewidth='4', linestyle='dotted', c='blue', label='SVC')
    ax3.minorticks_on()
    ax3.grid(which="both", alpha=0.3)
    ax3.legend()

    plt.tight_layout()
    plt.show()