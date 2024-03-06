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


epoch = 'july' # epoch to plot

lon = [296, 280] # degrees, lon positions to plot (always 2)
lat = [-21, -21] # degrees, lat positions to plot (always 2)

wave_min = [8.1, 7.2] # Range to confine plotted wavelengths to
wave_max = [11.0, 8.1]

fig_width = 12
fig_height = 9

shift = [2.0, 5.0] # Amount to shift the first coordinate data by

fig_name = 'figures/example_fit.png'


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

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(fig_width, fig_height))

colours = ['tab:blue', 'tab:red']
labels = ['a', 'b']

for kk in range(len_s):
	print('Plotting stage {}'.format(kk+1))

	hdul = fits.open(epoch + '/mosaics_stage{}/residual.fits'.format(kk+1))
	spec_cube = hdul['SPECTRAL_DATA'].data
	dspec_cube = hdul['ERR'].data
	fit_cube = hdul['FIT_DATA'].data
	wave_cube = hdul['WAVELENGTH'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	dz, dy, dx = spec_cube.shape

	wave = wave_cube[:, int(dy/2), int(dx/2)]

	min_idx = find_nearest(wave, wave_min[kk])
	max_idx = find_nearest(wave, wave_max[kk])

	wave = wave[min_idx:max_idx]
	spec_cube = spec_cube[min_idx:max_idx]
	dspec_cube = dspec_cube[min_idx:max_idx]
	fit_cube = fit_cube[min_idx:max_idx]

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

	for jj in range(len_p):
		print('Plotting position ({aa}, {bb})'.format(aa=lon[jj], bb=lat[jj]))

		lon_idx = find_nearest(lon_arr, lon[jj])
		lat_idx = find_nearest(lat_arr, lat[jj])

		if jj == 1:
			fit_cube += shift[kk]
			spec_cube += shift[kk]

		ax[kk].plot(wave, spec_cube[:, lat_idx, lon_idx], color=colours[jj], linestyle='-', linewidth=0.6, label='({aa}, {bb})'.format(aa=lon[jj], bb=lat[jj]))
		ax[kk].fill_between(wave, spec_cube[:, lat_idx, lon_idx]-dspec_cube[:, lat_idx, lon_idx], spec_cube[:, lat_idx, lon_idx]+dspec_cube[:, lat_idx, lon_idx], color=colours[jj], alpha=0.3)

		ax[kk].plot(wave, fit_cube[:, lat_idx, lon_idx], color='black', linestyle='-', linewidth=0.6)


	print('\n')

	ax[kk].set_xlabel('Wavelength (µm)')
	ax[kk].set_ylabel('Radiance (µWcm$^{-2}$sr$^{-1}$µm$^{-1}$)')
	ax[kk].set_title(r"$\bf{" + labels[kk] + "}$" + " Stage {}".format(kk+1), loc='left')

	ax[kk].grid()
	ax[kk].legend()

plt.subplots_adjust(hspace=0.3)

plt.savefig(fig_name, bbox_inches='tight')
plt.show()


#================================================================================

print('End of script\n')
