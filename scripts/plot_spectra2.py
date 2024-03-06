#================================================================================
# Description
"""
Plots comparison spectra for positions inputted below, generates the comparison
spectra for the GRS papers


Usage:
		[Input parameters below]

		python -W ignore plot_spectra.py


"""
#================================================================================
# Inputs


epoch = 'aug' # epoch to plot

lon = [310, 300] # degrees, lon positions to plot (always 2)
lat = [-16, -21] # degrees, lat positions to plot (always 2)

wave_min = 9.40 # Range to confine plotted wavelengths to
wave_max = 9.45

plot_ice = 9.42 # µm, wavelength to 2D plot

fig_width = 12
fig_height = 6

shift = 2.0 # Amount to shift the first coordinate data by

fig_name = 'figures/nh3_ice_search.png'


#================================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt

sys.path.insert(1, '/Users/jh852/Documents/project/python/subroutines')
from manifold import *


#================================================================================
# Main


len_s = 2
len_p = 2

print('\n')

# Data

hdul = fits.open(epoch + '/mosaics_stage1/residual.fits')
res_cube = hdul['RESIDUAL'].data
spec_cube = hdul['SPECTRAL_DATA'].data
dspec_cube = hdul['ERR'].data
fit_cube = hdul['FIT_DATA'].data
wave_cube = hdul['WAVELENGTH'].data
lon_grid = hdul['LON_WEST'].data
lat_grid = hdul['LAT_PGR'].data
hdul.close()

dz, dy, dx = spec_cube.shape

wave = wave_cube[:, int(dy/2), int(dx/2)]

min_idx = find_nearest(wave, wave_min)
max_idx = find_nearest(wave, wave_max)

wave = wave[min_idx:max_idx]
res_cube = res_cube[min_idx:max_idx]
spec_cube = spec_cube[min_idx:max_idx]
dspec_cube = dspec_cube[min_idx:max_idx]
fit_cube = fit_cube[min_idx:max_idx]

lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)



fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(fig_width, fig_height))


# Plot 2

xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

print(find_nearest(wave, plot_ice) + min_idx)
ax[1].imshow(res_cube[find_nearest(wave, plot_ice)], cmap='RdBu_r', origin='lower', vmin=-0.5, vmax=0.5)


ax[1].set_xlabel('West Longitude ($\degree$)')
ax[1].set_ylabel('Planetographic Latitude ($\degree$)')
ax[1].set_xticks(xtick_array, x_label_array)
ax[1].set_yticks(ytick_array, y_label_array)



# Plot 1

colours = ['tab:blue', 'tab:red']


for jj in range(len_p):
	print('Plotting position ({aa}, {bb})'.format(aa=lon[jj], bb=lat[jj]))

	lon_idx = find_nearest(lon_arr, lon[jj])
	lat_idx = find_nearest(lat_arr, lat[jj])

	ax[1].scatter(lon_idx, lat_idx, s=10, color='green')

	if jj == 1:
		fit_cube += shift
		spec_cube += shift

	ax[0].plot(wave, spec_cube[:, lat_idx, lon_idx], color=colours[jj], linestyle='-', linewidth=0.6, label='({aa}, {bb})'.format(aa=lon[jj], bb=lat[jj]))
	ax[0].fill_between(wave, spec_cube[:, lat_idx, lon_idx]-dspec_cube[:, lat_idx, lon_idx], spec_cube[:, lat_idx, lon_idx]+dspec_cube[:, lat_idx, lon_idx], color=colours[jj], alpha=0.3)

	ax[0].plot(wave, fit_cube[:, lat_idx, lon_idx], color='black', linestyle='-', linewidth=0.6)

ax[0].set_xlabel('Wavelength (µm)')
ax[0].set_ylabel('Residual (µWcm$^{-2}$sr$^{-1}$µm$^{-1}$)')

ax[0].grid()
ax[0].legend()






plt.savefig(fig_name, bbox_inches='tight')
plt.show()


#================================================================================

print('End of script\n')
