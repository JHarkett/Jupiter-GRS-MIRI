#======================================================================================
# Description
"""
Generates the main stratos figure

Usage
		[Input parameters below]

		python plot_stratos.py


"""
#======================================================================================
# Inputs


files = ['july/mosaics_stage2/000_retrieved_cleaned.fits', 'aug/mosaics_stage2/000_retrieved_cleaned.fits']
epoch_labels = ['2022-07-28', '2022-08-15']

alt_targ = [2, 8, 30] # mbar, targetted altitude for line plots

fig_width = 10
fig_height = 4

lon_c = [295.639, 301.0] # Lon centre of the GRS
lat_c = [-20.7, -20.7] # Lat centre of the GRS

lon_width = 11.09979 # Degrees, width of the GRS
lat_width = 8.880474 # Degrees, height of the GRS

max_pressure = 50.0 # mbar, min pressure of the x-sections
min_pressure = 0.1 # mbar, max pressure of the x-sections

label_pressures = [10.0, 1.0, 0.1] # Pressures to label

min_temp = 124 # K, min allowed temp in a and b
max_temp = 172 # K, max allowed temp in a and b

min_lat = -28.0 # Degrees, min lat to plot
max_lat = -13.0 # Degrees, max lat to plot

ref_file = 'nemesis.prf'
ref_skip = 24


#======================================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#======================================================================================
# Main


len_e = len(epoch_labels)

pressure = np.loadtxt(ref_file, skiprows=ref_skip, usecols=(1)) # atm
pressure *= 1013.25 # mbar

idx_p1 = find_nearest(pressure, max_pressure)
idx_p2 = find_nearest(pressure, min_pressure) + 1

pressure = pressure[idx_p1:idx_p2]

idx_targ1 = find_nearest(pressure, alt_targ[0])
idx_targ2 = find_nearest(pressure, alt_targ[1])
idx_targ3 = find_nearest(pressure, alt_targ[2])


# Generate the pressure ticks
pressure_denote_array = np.concatenate([np.arange(50.0, 10.0, -10.0), np.arange(10.0, 1.0, -1.0), np.arange(1.0, 0.0, -0.1)])
pressure_denote_array = np.round(pressure_denote_array, decimals=1)

pressuretick_array = [None]*len(pressure_denote_array)
pressure_label_array = [None]*len(pressure_denote_array)

print('\n')


for kk in range(len(pressure_denote_array)):

	pressuretick_array[kk] = find_nearest(pressure, pressure_denote_array[kk])

	if pressure_denote_array[kk] in label_pressures:
		print('Exception Triggered')

		pressure_label_array[kk] = pressure_denote_array[kk]
	else:
		pressure_label_array[kk] = None



# Plot configuration
fig = plt.figure(tight_layout=True, figsize=(fig_width, fig_height))
gs = gridspec.GridSpec(1, 3)

labels = ['a', 'b', 'c']

midi = (max_temp - min_temp)/2 + min_temp


for kk in range(len_e):
	print(epoch_labels[kk])

	hdul = fits.open(files[kk])
	temp_cube = hdul['VARIABLE'].data
	dtemp_cube = hdul['ERR'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

	lat_idx1 = find_nearest(lat_arr, min_lat)
	lat_idx2 = find_nearest(lat_arr, max_lat)

	lat_arr = lat_arr[lat_idx1:lat_idx2]
	xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

	zonal_temp = np.nanmedian(temp_cube, axis=2)
	zonal_dtemp = np.nanmedian(dtemp_cube, axis=2)

	zonal_temp = zonal_temp[idx_p1:idx_p2, lat_idx1:lat_idx2]
	zonal_dtemp = zonal_dtemp[idx_p1:idx_p2, lat_idx1:lat_idx2]

	dtemp_med = np.nanmedian(zonal_dtemp)


	# The zonal plot

	ax = fig.add_subplot(gs[kk])

	ax.imshow(zonal_temp, cmap='RdBu_r', origin='lower', vmin=min_temp, vmax=max_temp)

	ax.set_xlabel('Latitude ($\degree$)')
	ax.set_ylabel('Pressure (mbar)')
	ax.set_xticks(ytick_array, y_label_array)
	ax.set_yticks(pressuretick_array, pressure_label_array)

	ax.set_title(r"$\bf{" + labels[kk] + "}$" + " " + epoch_labels[kk], loc='left')

	# Contour lines
	ny, nx = zonal_temp.shape
	x, y = np.meshgrid(range(nx), range(ny))
	cplat = ax.contour(x, y, zonal_temp, np.arange(100, 200, 2), colors='black', linewidths=0.4, linestyles='solid')
	ax.clabel(cplat, fontsize='xx-small')


	# Colourbar
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = plt.colorbar(mappable=ax.imshow(zonal_temp, cmap='RdBu_r', origin='lower', vmin=min_temp, vmax=max_temp), cax=cax)
	cax.plot([0.5, 0.5], [midi-dtemp_med, midi+dtemp_med], 'k-', linewidth=1.0)
	cbar.ax.set_title('± {aa:.1f}'.format(aa=dtemp_med) + ' K', fontsize='x-small')


	# The line plot, note only plotting this for July, Aug is exactly the same


	if kk == 0:
		ax = fig.add_subplot(gs[2])

		# 1rd altitude
		ax.plot(lat_arr, zonal_temp[idx_targ1], color='tab:red', linestyle='-', linewidth=1.0, label=' {} mbar'.format(alt_targ[0]))
		ax.errorbar(lat_arr, zonal_temp[idx_targ1], yerr=zonal_dtemp[idx_targ1], xerr=[0.5]*len(lat_arr), fmt='none', ecolor='tab:red', elinewidth=0.7, capsize=3.0, capthick=0.7)

		# 2st altitude
		ax.plot(lat_arr, zonal_temp[idx_targ2], color='black', linestyle='-', linewidth=1.0, label=' {} mbar'.format(alt_targ[1]))
		ax.errorbar(lat_arr, zonal_temp[idx_targ2], yerr=zonal_dtemp[idx_targ2], xerr=[0.5]*len(lat_arr), fmt='none', ecolor='black', elinewidth=0.7, capsize=3.0, capthick=0.7)

		# 3nd altitude
		ax.plot(lat_arr, zonal_temp[idx_targ3], color='tab:blue', linestyle='-', linewidth=1.0, label=' {} mbar'.format(alt_targ[2]))
		ax.errorbar(lat_arr, zonal_temp[idx_targ3], yerr=zonal_dtemp[idx_targ2], xerr=[0.5]*len(lat_arr), fmt='none', ecolor='tab:blue', elinewidth=0.7, capsize=3.0, capthick=0.7)

		ax.set_xlabel('Latitude ($\degree$)')
		ax.set_ylabel('Retrieved Temperature (K)')
		ax.set_title(r"$\bf{" + 'c' + "}$" + " " + epoch_labels[0] + " Profile", loc='left')

		ax.grid()
		ax.legend(prop={'size':8})


plt.subplots_adjust(wspace=0.08)

plt.savefig('figures/stratos1_final.png', bbox_inches='tight')
plt.show()



#======================================================================================

print('End of script\n')
