#======================================================================================
# Description
"""
Generates the main stratos2 figure

Usage
		[Input parameters below]

		python plot_stratos2.py


"""
#======================================================================================
# Inputs


files1 = ['july/mosaics_stage2/000_retrieved_cleaned.fits', 'aug/mosaics_stage2/000_retrieved_cleaned.fits']
files2 = ['july/mosaics_stage2_aux/000_retrieved_cleaned.fits', 'aug/mosaics_stage2_aux/000_retrieved_cleaned.fits']

epoch_labels = ['2022-07-28', '2022-08-15']
stage_labels = ['Stage 2', 'Stage 2 AUX']

alt_targ = 3 # mbar, target altitude

fig_width = 10
fig_height = 5.5

lon_c = [295.639, 301.0] # Lon centre of the GRS
lat_c = [-20.7, -20.7] # Lat centre of the GRS

lon_extent = [37, 44] # Pixels, total width of FOV centred on the GRS
lat_extent = [19, 18] # Pixels, total height of FOV centred on the GRS

#lon_extent = [37, 1]   # Use these coordinates for
#lat_extent = [19, 18]  # analysing TMA 1

lon_width = 11.09979 # Degrees, width of the GRS
lat_width = 8.880474 # Degrees, height of the GRS

min_temp = 158 # K, min allowed temp
max_temp = 173 # K, max allowed temp


ref_file = 'nemesis.prf'
ref_skip = 24

analyse = False # If True: Will plot positions of the two TMA regions + give stats


#======================================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#======================================================================================
# Main


len_e = len(epoch_labels)
len_s = len(stage_labels)

pressure = np.loadtxt(ref_file, skiprows=ref_skip, usecols=(1)) # atm
pressure *= 1013.25 # mbar

idx_p = find_nearest(pressure, alt_targ)

print('\n')

# Plot configuration
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(fig_width, fig_height))

labels = ['a', 'b', 'c', 'd']

midi = (max_temp - min_temp)/2 + min_temp

for jj in range(2):
	print('Processing ' + stage_labels[jj] + '\n')

	if jj == 0:
		files = files1
	if jj == 1:
		files = files2

	for kk in range(len_e):
		print(epoch_labels[kk])

		hdul = fits.open(files[kk])
		temp_cube = hdul['VARIABLE'].data
		dtemp_cube = hdul['ERR'].data
		lon_grid = hdul['LON_WEST'].data
		lat_grid = hdul['LAT_PGR'].data
		hdul.close()

		lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

		lon_idx = find_nearest(lon_arr, lon_c[kk])
		lat_idx = find_nearest(lat_arr, lat_c[kk])

		lon_idx1 = lon_idx - lon_extent[0]
		lon_idx2 = lon_idx + lon_extent[1]
		lat_idx1 = lat_idx - lat_extent[0]
		lat_idx2 = lat_idx + lat_extent[1]

		temp_frame = temp_cube[idx_p, lat_idx1:lat_idx2, lon_idx1:lon_idx2]
		dtemp_frame = dtemp_cube[idx_p, lat_idx1:lat_idx2, lon_idx1:lon_idx2]
		lon_grid = lon_grid[lat_idx1:lat_idx2, lon_idx1:lon_idx2]
		lat_grid = lat_grid[lat_idx1:lat_idx2, lon_idx1:lon_idx2]

		dtemp_med = np.nanmedian(dtemp_frame)

		lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
		xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

		lon_idx = find_nearest(lon_arr, lon_c[kk])
		lat_idx = find_nearest(lat_arr, lat_c[kk])



		# The plot

		ax[jj, kk].imshow(temp_frame, cmap='inferno', origin='lower', vmin=min_temp, vmax=max_temp)

		ax[jj, kk].set_xlabel('West Longitude ($\degree$)')
		ax[jj, kk].set_ylabel('Latitude ($\degree$)')
		ax[jj, kk].set_xticks(xtick_array, x_label_array)
		ax[jj, kk].set_yticks(ytick_array, y_label_array)

		if jj == 0:
			ax[jj, kk].set_title(r"$\bf{" + labels[kk] + "}$" + " " + epoch_labels[kk] + " - " + stage_labels[jj] , loc='left')
		if jj == 1:
			ax[jj, kk].set_title(r"$\bf{" + labels[kk+2] + "}$" + " " + epoch_labels[kk] + " - " + stage_labels[jj] , loc='left')


		# Contour lines
		ny, nx = temp_frame.shape
		x, y = np.meshgrid(range(nx), range(ny))

		cplat = ax[jj, kk].contour(x, y, temp_frame, np.arange(120, midi, 1.5), colors='white', linewidths=0.4, linestyles='solid')
		ax[jj, kk].clabel(cplat, fontsize='xx-small')
		cplat = ax[jj, kk].contour(x, y, temp_frame, np.arange(midi, 190, 1.5), colors='black', linewidths=0.4, linestyles='solid')
		ax[jj, kk].clabel(cplat, fontsize='xx-small')


		# Colourbar
		divider = make_axes_locatable(ax[jj, kk])
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = plt.colorbar(mappable=ax[jj, kk].imshow(temp_frame, cmap='inferno', origin='lower', vmin=min_temp, vmax=max_temp), cax=cax)
		cax.plot([0.5, 0.5], [midi-dtemp_med, midi+dtemp_med], 'w-', linewidth=1.0)
		cbar.ax.set_title('± {aa:.1f}'.format(aa=dtemp_med) + ' K', fontsize='x-small')

		# GRS position boundaries
		elip = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_width*2, linewidth=2.0, linestyle='--', edgecolor='black', facecolor='none')
		ax[jj, kk].add_patch(elip)

		# TMA analysis

		if analyse == True and jj == 1:

			temp_tma = np.nanmax(temp_frame[y, x])

			print('Max Temp = {} K\n'.format(temp_tma))

			arral = np.where(temp_frame == temp_tma)

			print('Position = ({aa}, {bb})'.format(aa=lon_grid[arral[0][0], arral[1][0]], bb=lat_grid[arral[0][0], arral[1][0]]))
			print('Temp = {aa} ± {bb} K'.format(aa=temp_frame[arral[0][0], arral[1][0]], bb=dtemp_frame[arral[0][0], arral[1][0]]))
			print('Mean = {aa} ± {bb} K\n'.format(aa=np.nanmedian(temp_frame), bb=np.nanmedian(dtemp_frame)))

			ax[jj, kk].scatter(arral[1][0], arral[0][0], s=5, c='green')


plt.subplots_adjust(hspace=0.10, wspace=0.30)
plt.savefig('figures/stratos2_final.png', bbox_inches='tight')
plt.show()



#======================================================================================

print('End of script\n')
