#================================================================================
# Description
"""
Plots the Nyquist-combined data to check for any residual striping

Usage

		[Input parameters below]

		python plot_nyquist.py


"""
#================================================================================
# Inputs


file1 = 'aug/centre/d_combined/Level3_ch2-short_s3d_nav.fits' # Original
file2 = 'aug/centre/d_combined_0.5/Level3_ch2-short_s3d_nav.fits' # New

idx_plot = 1 # Starting from 1
wave_plot = 7.51065 # µm

lat_targ = -24.0 # Latitude

fig_width = 8
fig_height = 8


#================================================================================
# Imports


import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#================================================================================
# Main


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_width, fig_height))
colours = ['black', 'red']
labels = ['Original (stripy)', 'New']


for kk in range(2):
	if kk == 0:
		hdul = fits.open(file1)

	if kk == 1:
		hdul = fits.open(file2)

	spec_cube = hdul['SCI'].data
	dz, dy, dx = spec_cube.shape

	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

	y_pos = find_nearest(lat_arr, lat_targ)

	hdul.close()

	spec = spec_cube[idx_plot-1, y_pos, :]

	lat_val = lat_arr[y_pos]

	ax.plot(lon_arr, spec, linestyle='-', linewidth=0.8, color=colours[kk], label=labels[kk] + ' lat={}'.format(lat_val))


ax.set_xlabel('West Longitude ($\degree$)')
ax.set_ylabel('Brightness (MJy/sr)')

ax.grid()
ax.set_xlim(max(lon_arr), min(lon_arr))
ax.legend()
ax.set_title('Wavelength = {aa} µm, idx = {bb}, y = {cc}'.format(aa=wave_plot, bb=idx_plot, cc=lat_targ))


plt.show()


#================================================================================

print('End of script\n')
