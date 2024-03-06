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


files = ['july/mosaics_stage1/000_retrieved_cleaned.fits', 'aug/mosaics_stage1/000_retrieved_cleaned.fits']

epoch_labels = ['2022-07-28', '2022-08-15']

alt_targ = [100, 300, 400, 500, 700] # mbar
units = ['K', 'K', 'K', 'K', 'K'] # Units on the colourbar

min_allow = [98, 110, 117, 128, 141] # min allowed value
max_allow = [114, 123, 128, 139, 149] # max allowed value

intervals = [3, 2, 2, 2, 2] # Interval of the contour lines

fig_width = 8
fig_height = 9

fig_name = 'figures/tropos_final.png'

trim_july = [[31, 2], [2, 0]] # [[x1, x2], [y1, y2]]
trim_aug = [[1, 33], [3, 1]] # [[x1, x2], [y1, y2]]

lon_c = [295.639, 301.0] # Lon centre of the GRS
lat_c = [-20.7, -20.7] # Lat centre of the GRS

lon_width = 11.09979 # Degrees, width of the GRS
lat_width = 8.880474 # Degrees, height of the GRS

lon_storm = 309.5 # W lon of the GRS wake storm
lat_storm = -15.0 # PGR lat of said storm


#=====================================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=====================================================================================
# Main


len_e = len(epoch_labels)
len_s = len(alt_targ) # Number of altitudes to plot

fig = plt.figure(figsize=(fig_width, fig_height))
gs = gridspec.GridSpec(5, 4, width_ratios=(2, 2.2, 1, 2))

print('\n')

data_pre = np.loadtxt('tempapr.dat', skiprows=1, usecols=(0, 1))

pressure = data_pre[:,0] # atm
pressure *= 1013.25 # mbar

temp_prior = data_pre[:,1] # K

colours22000 = ['black', 'tab:red']
ax22000 = fig.add_subplot(gs[:, 3])

idx_deep = find_nearest(pressure, 1000)
idx_shallow = find_nearest(pressure, 50)
pressure2 = pressure[idx_deep:idx_shallow]

ax22000.plot([0, 0], [1000, 50], color='gray', linestyle='--', linewidth=1.0)


for kk in range(len_e):
	print('Opening ' + epoch_labels[kk] + ' data\n')

	# Temp data

	hdul = fits.open(files[kk])
	temp = hdul['VARIABLE'].data
	dtemp = hdul['ERR'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	dz, dy, dx = temp.shape

	if kk == 0:
		temp = temp[:, trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		dtemp = dtemp[:, trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lon_grid = lon_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lat_grid = lat_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]

	if kk == 1:
		temp = temp[:, trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		dtemp = dtemp[:, trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lon_grid = lon_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lat_grid = lat_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
	xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

	lon_idx = find_nearest(lon_arr, lon_c[kk])
	lat_idx = find_nearest(lat_arr, lat_c[kk])



	# Plot 2

	x_int = np.arange(0, dx)
	y_int = np.arange(0, dy)

	oval = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_width*2, linewidth=2.0, linestyle='--', edgecolor='tab:green', facecolor='none')

	g = np.meshgrid(x_int, y_int)
	coords = list(zip(*(c.flat for c in g)))

	ellipsepoints_in = np.vstack([p for p in coords if oval.contains_point(p, radius=0)])

	in_x = ellipsepoints_in[:,0]
	in_y = ellipsepoints_in[:,1]

	temp_mean_grs = np.nanmean(temp[:, in_y, in_x], axis=1)
	err_mean_grs = np.nanmean(dtemp[:, in_y, in_x], axis=1)
	dtemp_mean_grs = temp_mean_grs - temp_prior

	dtemp_mean_grs = dtemp_mean_grs[idx_deep:idx_shallow]
	err_mean_grs = err_mean_grs[idx_deep:idx_shallow]


	ax22000.plot(dtemp_mean_grs, pressure2, color=colours22000[kk], linestyle='-', linewidth=0.8, label=epoch_labels[kk])
	ax22000.fill_betweenx(pressure2, dtemp_mean_grs-err_mean_grs, dtemp_mean_grs+err_mean_grs, color=colours22000[kk], alpha=0.3)


	if kk == 1:
		idx_x_storm = find_nearest(lon_arr, lon_storm)
		idx_y_storm = find_nearest(lat_arr, lat_storm)

		oval_storm = patches.Ellipse((idx_x_storm, idx_y_storm), width=8, height=8, linewidth=2.0, linestyle='--', edgecolor='tab:green', facecolor='none')

		ellipsepoints_storm = np.vstack([p for p in coords if oval_storm.contains_point(p, radius=0)])
		storm_x = ellipsepoints_storm[:,0]
		storm_y = ellipsepoints_storm[:,1]

		temp_mean_storm = np.nanmean(temp[:, storm_y, storm_x], axis=1)
		err_mean_storm = np.nanmean(dtemp[:, storm_y, storm_x], axis=1)
		dtemp_mean_storm = temp_mean_storm - temp_prior

		dtemp_mean_storm = dtemp_mean_storm[idx_deep:idx_shallow]
		err_mean_storm = err_mean_storm[idx_deep:idx_shallow]

		ax22000.plot(dtemp_mean_storm, pressure2, color='tab:blue', linestyle='-', linewidth=0.8, label='VICI 1')
		ax22000.fill_betweenx(pressure2, dtemp_mean_storm-err_mean_storm, dtemp_mean_storm+err_mean_storm, color='blue', alpha=0.3)



	# The original plot
	for jj in range(len_s):
		idx_use = find_nearest(pressure, alt_targ[jj])

		plot_species = temp[idx_use]
		plot_err = dtemp[idx_use]

		ax = fig.add_subplot(gs[jj, kk])

		ax.imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj])

		elip = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_width*2, linewidth=2.0, linestyle='--', edgecolor='tab:green', facecolor='none')
		ax.add_patch(elip)

		ax.set_xticks(xtick_array, [])
		ax.set_yticks(ytick_array, [])

		if kk == 0:
			ax.set_yticks(ytick_array, y_label_array)
			ax.tick_params(labelsize='small')

		if jj == 4:
			ax.set_xticks(xtick_array, x_label_array)
			ax.tick_params(labelsize='small')


		# Contour lines

		ny, nx = plot_species.shape
		x, y = np.meshgrid(range(nx), range(ny))

		if min_allow[jj] == None:
			midi = (np.nanmax(plot_species) - np.nanmin(plot_species))/2 + np.nanmin(plot_species)

			cplat = ax.contour(x, y, plot_species, np.arange(midi, np.nanmax(plot_species), intervals[jj]), colors='black', linewidths=0.4, linestyles='solid')
			ax.clabel(cplat, fontsize='xx-small')

			cplat = ax.contour(x, y, plot_species, np.arange(np.nanmin(plot_species), midi, intervals[jj]), colors='white', linewidths=0.4, linestyles='solid')
			ax.clabel(cplat, fontsize='xx-small')

		else:
			midi = (max_allow[jj] - min_allow[jj])/2 + min_allow[jj]

			cplat = ax.contour(x, y, plot_species, np.arange(midi, max_allow[jj], intervals[jj]), colors='black', linewidths=0.4, linestyles='solid')
			ax.clabel(cplat, fontsize='xx-small')

			cplat = ax.contour(x, y, plot_species, np.arange(min_allow[jj], midi, intervals[jj]), colors='white', linewidths=0.4, linestyles='solid')
			ax.clabel(cplat, fontsize='xx-small')

		if kk == 1:

			err_median = np.nanmedian(plot_err)

			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			cbar = plt.colorbar(mappable=ax.imshow(plot_species, cmap='inferno', origin='lower', vmin=min_allow[jj], vmax=max_allow[jj]), cax=cax)
			cbar.ax.set_title('± {aa:.2f} '.format(aa=err_median) + units[jj], fontsize='x-small')

			if min_allow[jj] != None:
				midi = (max_allow[jj] - min_allow[jj])/2 + min_allow[jj]
				plt.plot([0.5, 0.5], [midi-err_median, midi+err_median], 'w-', linewidth=1.0)

			oval_storm = patches.Ellipse((idx_x_storm, idx_y_storm), width=8, height=8, linewidth=2.0, linestyle='--', edgecolor='tab:cyan', facecolor='none')
			ax.add_patch(oval_storm)

		if kk == 0:
			ax.set_ylabel('{} mbar'.format(alt_targ[jj]))

		if jj == 0:
			ax.set_title(epoch_labels[kk])

ax22000.set_title('a                                                                              b', weight='bold', loc='left', x=-2.85, y=1.01)
ax22000.legend()
ax22000.grid()
ax22000.set_xlabel('$\Delta$T (K)')
ax22000.set_ylabel('Pressure (mbar)')
ax22000.set_ylim(max(pressure2), min(pressure2))


fig.supxlabel('West Longitude ($\degree$)', x=0.35, y=0.05, fontsize='medium')
fig.supylabel('Planetograpic Latitude ($\degree$)', x=0.03, fontsize='medium')

plt.subplots_adjust(wspace=0.1, hspace=0.0)

plt.savefig(fig_name, bbox_inches='tight')
plt.show()


#=====================================================================================

print('End of script\n')
