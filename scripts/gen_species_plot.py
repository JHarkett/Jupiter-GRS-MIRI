#=====================================================================================
# Description
"""
Loads and plots the stratospheric altitudes input below

Usage

		[Input parameters below]

		python gen_species_plot.py


"""
#=====================================================================================
# Inputs


nh3_files = ['july/mosaics_stage1/11020_retrieved_cleaned.fits', 'aug/mosaics_stage1/11020_retrieved_cleaned.fits']
ph3_files = ['july/mosaics_stage1/28020_retrieved_cleaned.fits', 'aug/mosaics_stage1/28020_retrieved_cleaned.fits']
aero_files = ['july/mosaics_stage1/-1032_retrieved_cleaned.fits', 'aug/mosaics_stage1/-1032_retrieved_cleaned.fits']

epoch_labels = ['2022-07-28', '2022-08-15']

row_labels = ['NH$_3$ VMR', 'NH$_3$ FSH', 'PH$_3$ VMR', 'PH$_3$ FSH', 'Aerosol Opacity']
units = ['ppm', '', 'ppm', '', ''] # Units on the colourbar

min_allow = [40, 0.08, 0.3, 0.03, 0.9] # min allowed value
max_allow = [450, 0.15, 1.0, 0.2, 4.0] # max allowed value

intervals = [75, 0.01, 0.1, 0.04, 0.6] # Interval of the contour lines

nama = ['Ammonia at 800 mbar', 'Phosphine at 500 mbar', 'Aerosol Opacity']


fig_width = 6
fig_height = 8

fig_name = 'figures/species_final.png'

trim_july = [[31, 2], [2, -1]] # [[x1, x2], [y1, y2]]
trim_aug = [[1, 33], [3, 2]] # [[x1, x2], [y1, y2]]

lon_c = [295.639, 301.0] # Lon centre of the GRS
lat_c = [-20.7, -20.7] # Lat centre of the GRS

lon_width = 11.09979 # Degrees, width of the GRS
lat_width = 8.880474 # Degrees, height of the GRS


#=====================================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=====================================================================================
# Main


len_e = len(epoch_labels)
len_s = 5 # Number of species parameters to plot

fig, ax = plt.subplots(nrows=len_s, ncols=len_e, figsize=(fig_width, fig_height), tight_layout=True, sharey=True, sharex='col')

print('\n')


# Aux plot 2
#fig2, ax2 = plt.subplots(nrows=3, ncols=1, figsize=(4, 7), squeeze=False)



for kk in range(len_e):
	print('Opening ' + epoch_labels[kk] + ' data\n')

	# NH3 data

	hdul = fits.open(nh3_files[kk])
	nh3_vmr = hdul['VAL1'].data * 1e+6
	nh3_dvmr = hdul['ERR1'].data * 1e+6
	nh3_fsh = hdul['VAL2'].data
	nh3_dfsh = hdul['ERR2'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	# PH3 data

	hdul = fits.open(ph3_files[kk])
	ph3_vmr = hdul['VAL1'].data * 1e+6
	ph3_dvmr = hdul['ERR1'].data * 1e+6
	ph3_fsh = hdul['VAL2'].data
	ph3_dfsh = hdul['ERR2'].data
	hdul.close()

	# Aerosol data

	hdul = fits.open(aero_files[kk])
	aero_op = hdul['VAL2'].data
	aero_dop = hdul['ERR2'].data
	hdul.close()

	dy, dx = nh3_vmr.shape


	if kk == 0:
		nh3_vmr = nh3_vmr[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		ph3_vmr = ph3_vmr[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		nh3_fsh = nh3_fsh[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		ph3_fsh = ph3_fsh[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		aero_op = aero_op[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lon_grid = lon_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lat_grid = lat_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]

	if kk == 1:
		nh3_vmr = nh3_vmr[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		ph3_vmr = ph3_vmr[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		nh3_fsh = nh3_fsh[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		ph3_fsh = ph3_fsh[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		aero_op = aero_op[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lon_grid = lon_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lat_grid = lat_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]


	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
	xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

	lon_idx = find_nearest(lon_arr, lon_c[kk])
	lat_idx = find_nearest(lat_arr, lat_c[kk])


	# The plot
	for jj in range(len_s):
		if jj == 0:
			plot_species = nh3_vmr
			plot_err = nh3_dvmr
		if jj == 1:
			plot_species = nh3_fsh
			plot_err = nh3_dfsh
		if jj == 2:
                        plot_species = ph3_vmr
                        plot_err = ph3_dvmr
		if jj == 3:
			plot_species = ph3_fsh
			plot_err = ph3_dfsh
		if jj == 4:
			plot_species = aero_op
			plot_err = aero_dop

		ax[jj, kk].imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj])

		print(np.nanmean(plot_species))

		elip = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_width*2, linewidth=2.0, linestyle='--', edgecolor='red', facecolor='none')
		ax[jj, kk].add_patch(elip)

		ax[jj, kk].set_xticks(xtick_array, x_label_array)
		ax[jj, kk].set_yticks(ytick_array, y_label_array)


		# Contour lines

		midi = (max_allow[jj] - min_allow[jj])/2 + min_allow[jj]

		ny, nx = plot_species.shape
		x, y = np.meshgrid(range(nx), range(ny))
		cplat = ax[jj, kk].contour(x, y, plot_species, np.arange(midi, max_allow[jj], intervals[jj]), colors='black', linewidths=0.4, linestyles='solid')
		ax[jj, kk].clabel(cplat, fontsize='xx-small')

		cplat = ax[jj, kk].contour(x, y, plot_species, np.arange(min_allow[jj], midi, intervals[jj]), colors='white', linewidths=0.4, linestyles='solid')
		ax[jj, kk].clabel(cplat, fontsize='xx-small')


		if kk == 1:

			err_median = np.nanmedian(plot_err)

			divider = make_axes_locatable(ax[jj, kk])
			cax = divider.append_axes("right", size="5%", pad=0.05)
			cbar = plt.colorbar(mappable=ax[jj, kk].imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj]), cax=cax)

			if jj == 1 or jj == 3:
				cbar.ax.set_title('± {aa:.3f} '.format(aa=err_median) + units[jj], fontsize='x-small')
			else:
				cbar.ax.set_title('± {aa:.2f} '.format(aa=err_median) + units[jj], fontsize='x-small')

			midi = (max_allow[jj] - min_allow[jj])/2 + min_allow[jj]

			plt.plot([0.5, 0.5], [midi-err_median, midi+err_median], 'w-', linewidth=1.0)


				# Aux plot 2

				#ax2[jj-int(jj/2), kk-1].imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj])

				#elip = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_width*2, linewidth=2.0, linestyle='--', edgecolor='red', facecolor='none')
				#ax2[jj-int(jj/2), kk-1].add_patch(elip)

				#ax2[jj-int(jj/2), kk-1].set_xticks(xtick_array, x_label_array)
				#ax2[jj-int(jj/2), kk-1].set_yticks(ytick_array, y_label_array)

				#ax2[jj-int(jj/2), kk-1].set_title(nama[jj-int(jj/2)])

				#err_median = np.nanmedian(plot_err)

				#divider = make_axes_locatable(ax2[jj-int(jj/2), kk-1])
				#cax = divider.append_axes("right", size="5%", pad=0.05)
				#cbar = plt.colorbar(mappable=ax2[jj-int(jj/2), kk-1].imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj]), cax=cax)
				#cbar.ax.set_title('± {aa:.2f} '.format(aa=err_median) + units[jj], fontsize='x-small')

				#midi = (max_allow[jj] - min_allow[jj])/2 + min_allow[jj]

				#plt.plot([0.5, 0.5], [midi-err_median, midi+err_median], 'k-', linewidth=1.0)




		if kk == 0:
			ax[jj, kk].set_ylabel(row_labels[jj])

		if jj == 0:
			ax[jj, kk].set_title(epoch_labels[kk])


#fig2.supxlabel('Longitude ($\degree$)', y=0.04)
#fig2.supylabel('Latitude ($\degree$)', x=0.05)
#plt.subplots_adjust(hspace=0.4)

#plt.savefig('figures/species_poster.png', bbox_inches='tight')
#plt.show()



fig.supxlabel('West Longitude ($\degree$)', y=0.02)
fig.supylabel('Planetograpic Latitude ($\degree$)', x=0.02)

plt.subplots_adjust(wspace=-1.75, hspace=-1.1)

plt.savefig(fig_name, bbox_inches='tight')
plt.show()




#=====================================================================================

print('End of script\n')
