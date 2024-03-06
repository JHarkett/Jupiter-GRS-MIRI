#====================================================================================
# Description
"""
Plots coords on a FOV input below

Usage
		[Input parameters below]

		python plot_coords.py

"""
#====================================================================================
# Inputs


file = 'aug/west/d_combined_0.5/Level3_ch2-short_s3d_nav.fits' # File to plot
wave_plot = 8.6 # µm, wavelength to visulise

lon = 313.0 # Degrees, target longitude
lat = -15.0 # Degrees, target latitude

fig_width = 8
fig_height = 7


#====================================================================================
# Imports


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#====================================================================================
# Main


hdul = fits.open(file)
spec_cube = hdul['SCI'].data
lon_grid = hdul['LON_WEST'].data
lat_grid = hdul['LAT_PGR'].data
hdul.close()

wave = gen_wave_arr(file)

lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

wave_idx = find_nearest(wave, wave_plot)
spec_frame = spec_cube[wave_idx]

lon_idx = find_nearest(lon_arr, lon)
lat_idx = find_nearest(lat_arr, lat)


# The plot

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_width, fig_height))

ax.imshow(spec_frame, cmap='inferno', origin='lower')
ax.scatter(lon_idx, lat_idx, s=20, c='green')

ax.set_xlabel('West Longitude ($\degree$)')
ax.set_ylabel('PGR Latitude ($\degree$)')
ax.set_xticks(xtick_array, x_label_array)
ax.set_yticks(ytick_array, y_label_array)
ax.set_title('{aa} µm frame, ({bb}, {cc})'.format(aa=wave_plot, bb=lon, cc=lat))

plt.show()


#====================================================================================

print('End of script\n')
