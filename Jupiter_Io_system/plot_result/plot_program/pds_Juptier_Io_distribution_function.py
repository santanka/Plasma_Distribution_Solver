import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants


plt.rcParams.update({'mathtext.default': 'default', 'mathtext.fontset': 'stix'})
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 60

data = np.genfromtxt(r"/home/satanka/Documents/Plasma_Distribution_Solver/Jupiter_Io_system/result_distribution_function/result_distribution_function_027_027_BC_002_min_063_grid_050_series_06.csv", delimiter=',', unpack=True)

channel = 1

mlat_degree = data[0, :]
v_perp_i = data[1, :]
v_para_i = data[2, :]
distribution_function_i = data[3, :]
v_perp_b = data[4, :]
v_para_b = data[5, :]
distribution_function_b = data[6, :]

if (channel == 1):
    length = len(distribution_function_i)
    max_distribution_function_i = np.nanmax(distribution_function_i)
    for ii in range(length):
        if(np.log10(distribution_function_i[ii]) < np.log10(max_distribution_function_i)-20.):
            v_perp_i[ii] = np.nan
            v_para_i[ii] = np.nan
            distribution_function_i[ii] = np.nan
            v_perp_b[ii] = np.nan
            v_para_b[ii] = np.nan
            distribution_function_b[ii] = np.nan

    fig = plt.figure()
    ax = fig.add_subplot(111)

    cm = plt.cm.get_cmap(r'turbo')

    ax.set_xlabel(r"$v_{\parallel}$ $[\mathrm{m/s}]$ (+ : South → North, - : North → South)")
    ax.set_ylabel(r"$v_{\perp}$ $[\mathrm{m/s}]$")
    plt.title(r"Distribution Function (scale=log10)")

    if(min(np.floor(distribution_function_i)) != 0.):
        mappable = ax.scatter(v_para_i, v_perp_i, c=np.log10(distribution_function_i), vmin=np.floor(np.nanmin(np.log10(distribution_function_i))), vmax=np.trunc(np.nanmax(np.log10(distribution_function_i))), cmap=cm, s=700, alpha=0.7)
    if(min(np.floor(distribution_function_i)) == 0.):
        mappable = ax.scatter(v_para_i, v_perp_i, c=np.log10(distribution_function_i), vmin=np.floor(np.nanmax(np.log10(distribution_function_i))-15.), vmax=np.trunc(np.nanmax(np.log10(distribution_function_i))), cmap=cm, s=700, alpha=0.7)

    cbar = fig.colorbar(mappable, ax=ax, label=r"$\log_{10} f_{si}$")

if (channel == 2):
    length = len(distribution_function_b)
    max_distribution_function_b = np.nanmax(distribution_function_b)
    for ii in range(length):
        if(np.log10(distribution_function_b[ii]) < np.log10(max_distribution_function_b)-20.):
            v_perp_i[ii] = np.nan
            v_para_i[ii] = np.nan
            distribution_function_i[ii] = np.nan
            v_perp_b[ii] = np.nan
            v_para_b[ii] = np.nan
            distribution_function_b[ii] = np.nan

    fig = plt.figure()
    ax = fig.add_subplot(111)

    cm = plt.cm.get_cmap(r'turbo')

    ax.set_xlabel(r"$v_{\parallel}$ $[\mathrm{m/s}]$ (+ : South → North, - : North → South)")
    ax.set_ylabel(r"$v_{\perp}$ $[\mathrm{m/s}]$")
    plt.title(r"Distribution Function (scale=log10)")

    if(min(np.floor(distribution_function_b)) != 0.):
        mappable = ax.scatter(v_para_b, v_perp_b, c=np.log10(distribution_function_b), vmin=np.floor(np.nanmin(np.log10(distribution_function_b))), vmax=np.trunc(np.nanmax(np.log10(distribution_function_b))), cmap=cm, s=700, alpha=0.7)
    if(min(np.floor(distribution_function_i)) == 0.):
        mappable = ax.scatter(v_para_b, v_perp_b, c=np.log10(distribution_function_b), vmin=np.floor(np.nanmax(np.log10(distribution_function_b))-15.), vmax=np.trunc(np.nanmax(np.log10(distribution_function_b))), cmap=cm, s=700, alpha=0.7)

    cbar = fig.colorbar(mappable, ax=ax, label=r"$\log_{10} f_{sb}$")

ax.minorticks_on()
ax.grid(which="both")
ax.set_axisbelow(True)
plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.show()