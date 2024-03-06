#=====================================================================================
# Description
"""
Loads the temperature data from input files and pressures
below, plots 3D temperature maps

Usage:

		[Input parameters below]

		python plot_temperatures.py


"""
#=====================================================================================
# Inputs


# List of files that contain the temperature data
#	Temperature should be under a fits header called 'VARIABLE'
#	Uncertainty should be under a fits header called 'ERR'
#	Longitude data should be under a fits header called 'LON_WEST'
#	Latitude data should be under a fits header called 'LAT_PGR'

files = ['july/mosaics_stage2_aux/000_retrieved_cleaned.fits',\
'aug/mosaics_stage2_aux/000_retrieved_cleaned.fits']


titles = 'July', 'August'

pressure_file = 'jupiter_v2022miri.ref'  # Location of the pressure array (atm)
ref_skip = 25  # number of rows to skip to get to the pressure data

# mbar
pressures_use = [1.0, 3.0, 8.0, 30.0, 'HST']

pressures_names = ['1.0 mbar', '3.0 mbar', '8.0 mbar', '30.0 mbar', 'Context']

stretch = 80

elev = 25
azim = 230
roll = None


trim_july = [[30, 2], [0, 0]]
trim_aug = [[0, 10], [1, 1]]

# Location of context data
hst_file = 'context_data/hst_fig1.png'
alpo_file = 'context_data/alpo_aug.png'

save_name = 'figures/tropos_strat3d.png' # Name of file to create

fig_width = 10
fig_height = 5


#=====================================================================================
# Imports


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import scipy.interpolate as si

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=====================================================================================
# Functions


def NormalizeData(data):
	return (data - (np.nanmin(data)-1)) / ((np.nanmax(data)+1) - (np.nanmin(data)-1))


class Arrow3D(FancyArrowPatch):
	def __init__(self, xs, ys, zs, *args, **kwargs):
		super().__init__((0,0), (0,0), *args, **kwargs)
		self._verts3d = xs, ys, zs

	def do_3d_projection(self, renderer=None):
		xs3d, ys3d, zs3d = self._verts3d
		xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
		self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

		return np.min(zs)


#=====================================================================================
# Main



len_p = len(pressures_names)
len_f = len(files)

pressure = np.loadtxt(pressure_file, skiprows=ref_skip, usecols=(1))  # atm
pressure *= 1013.25  # mbar

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(fig_width, fig_height), subplot_kw={'projection':'3d'})

markers = ['a', 'b', 'c', 'd']

for jj in range(len_f):

	hdul = fits.open(files[jj])
	temp = hdul['VARIABLE'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	dz, dy, dx = temp.shape

	if jj == 0:
		temp = temp[:, trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lon_grid = lon_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]
		lat_grid = lat_grid[trim_july[1][0]:dy-trim_july[1][1], trim_july[0][0]:dx-trim_july[0][1]]

	if jj == 1:
		temp = temp[:, trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lon_grid = lon_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]
		lat_grid = lat_grid[trim_aug[1][0]:dy-trim_aug[1][1], trim_aug[0][0]:dx-trim_aug[0][1]]

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)


	dz, dy, dx = temp.shape


	if jj == 0:
		hst, lon_arr_hst, lat_arr_hst, lon_grid_hst, lat_grid_hst = open_png(hst_file, 315.0, 274.0, -30.0, -10.0)
	if jj == 1:
		hst, lon_arr_hst, lat_arr_hst, lon_grid_hst, lat_grid_hst = open_png(alpo_file, 369.0, 158.0, -45.0, -5.0)


	lat1 = find_nearest(lat_arr_hst, np.nanmin(lat_arr))
	lat2 = find_nearest(lat_arr_hst, np.nanmax(lat_arr))
	lon1 = find_nearest(lon_arr_hst, np.nanmax(lon_arr))
	lon2 = find_nearest(lon_arr_hst, np.nanmin(lon_arr))

	hst = hst[lat1:lat2, lon1:lon2, :]
	lon_grid_hst = lon_grid_hst[lat1:lat2, lon1:lon2]
	lat_grid_hst = lat_grid_hst[lat1:lat2, lon1:lon2]
	lon_arr_hst = lon_arr_hst[lon1:lon2]
	lat_arr_hst = lat_arr_hst[lat1:lat2]


	B_frame = hst[:, :, 0]
	G_frame = hst[:, :, 1]
	R_frame = hst[:, :, 2]

	B_frame = NormalizeData(B_frame)
	G_frame = NormalizeData(G_frame)
	R_frame = NormalizeData(R_frame)

	dy1, dx1, dz1 = hst.shape

	print(dx1, dy1, dz1)
	print(dx, dy, dz)

	rgbArray = np.zeros((dy1, dx1, 3))
	rgbArray[:, :, 0] = B_frame
	rgbArray[:, :, 1] = G_frame
	rgbArray[:, :, 2] = R_frame


	x_grid, y_grid = np.meshgrid(np.arange(0, dx), np.arange(0, dy))

	z_poss = [60, 50, 40, 30, 20, 10, 0]

	x_grid_hst, y_grid_hst = np.meshgrid(np.arange(0, dx1), np.arange(0, dy1))

	min_temps = []
	max_temps = []

	for kk in range(len_p):
		print(kk)

		z_grid = np.ones((dy, dx)) * z_poss[kk]

		if pressures_use[kk] == 'HST':
			if jj == 0:
				ax[jj].plot_surface(x_grid, y_grid, z_grid, facecolors = rgbArray[::20, ::20, :],\
				antialiased = True, rstride=1, cstride=1, alpha=None)
			else:
				ax[jj].plot_surface(x_grid, y_grid, z_grid, facecolors = rgbArray[::4, ::4, :],\
				antialiased = True, rstride=1, cstride=1, alpha=None)

		else:
			temp_layer = temp[find_nearest(pressure, pressures_use[kk])]

			scam = plt.cm.ScalarMappable(norm=cm.colors.Normalize(np.nanmin(temp_layer), np.nanmax(temp_layer)), cmap='inferno')

			ax[jj].plot_surface(x_grid, y_grid, z_grid, facecolors = scam.to_rgba(temp_layer),\
			antialiased = True, rstride=1, cstride=1, alpha=None)

		ax[jj].text(dx, 0, z_poss[kk], pressures_names[kk])




	a = Arrow3D([dx+5, dx+5], [15, 35], [z_poss[0], z_poss[0]], mutation_scale=10, lw=2, arrowstyle="-|>", color="k")
	ax[jj].add_artist(a)
	ax[jj].text(dx+10, 25, z_poss[0], 'North')

	ax[jj].view_init(elev=elev, azim=azim, roll=roll, vertical_axis='z')
	ax[jj].set_box_aspect([dx, dy, stretch])
	ax[jj].set_axis_off()

	ax[jj].set_title(r"$\bf{" + markers[jj] + "}$" + " " + titles[jj])


cbar_ax = fig.add_axes([0.4, 0.19, 0.2, 0.02])
fig.colorbar(mappable=scam, cax=cbar_ax, ticks=[], orientation='horizontal', label='Cooler        Warmer')
cbar_ax.arrow(120.0, 0.5, -4.5, 0.0, length_includes_head=True, head_width=0.4, color='white')
cbar_ax.arrow(129.0, 0.5, 4.5, 0.0, length_includes_head=True, head_width=0.4, color='black')
cbar_ax.tick_params(labelsize=5)


plt.subplots_adjust(wspace=0.00)
plt.savefig(save_name, bbox_inches='tight')
plt.show()



#=====================================================================================

print('End of script\n')


