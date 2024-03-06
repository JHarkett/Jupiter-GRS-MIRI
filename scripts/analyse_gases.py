#===========================================================================
# Description
"""
Compares NH3 FSH to expected SVP
May also eventually also compare PH3 deep VMR to aerosol opacity too


Usage:
		[Input parameters below]

		python analyse_gases.py


"""
#===========================================================================
# Inputs


epochs = ['july', 'aug']

nh3_file = 'mosaics_stage1/11020_retrieved_cleaned.fits'
ph3_file = 'mosaics_stage1/28020_retrieved_cleaned.fits'
aero_file = 'mosaics_stage1/-1032_retrieved_cleaned.fits'
temp_file = 'mosaics_stage1/000_retrieved_cleaned.fits'

nh3_alt = 800 # Assumed altitude of the NH3 cloud layer (mbar)

fig_width = 11
fig_height = 10

lon_c = [295.639, 301.0] # Lon centre of the GRS
lat_c = [-20.7, -20.7] # Lat centre of the GRS

lon_width = 11.09979
lat_height = 8.880474

pressure_file = 'jupiter_v2022miri.ref'
ref_skip = 25

wspace = 0.2
hspace = 0.3

fig_name = 'figures/gaseous_analysis_recent_final.png'


#===========================================================================
# Imports and useful values


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *

# NH3 SVP constants
A = 23.224
B = -4245.8
C = -2.2775E-2
# Note: D == 0.0 in this case so is not included in the calculation


#===========================================================================
# Main


pressure = np.loadtxt(pressure_file, skiprows=ref_skip, usecols=(1)) # atm
pressure *= 1013.25 # mbar


fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(fig_width, fig_height))

labels = [['a', 'b'], ['c', 'd']]
colours = ['black', 'tab:red']
sub_labels = ['GRS', 'Outside']

for cccp in range(2):

	hdul = fits.open(epochs[cccp] + '/' + temp_file)
	temp = hdul['VARIABLE'].data
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	dz, dy, dx = temp.shape

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
	xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

	lon_idx = find_nearest(lon_arr, lon_c[cccp])
	lat_idx = find_nearest(lat_arr, lat_c[cccp])


	oval = patches.Ellipse((lon_idx, lat_idx), lon_width*2, lat_height*2)

	x_int = np.arange(0, dx)
	y_int = np.arange(0, dy)

	g = np.meshgrid(x_int, y_int)
	coords = list(zip(*(c.flat for c in g)))

	ellipsepoints_in = np.vstack([p for p in coords if oval.contains_point(p, radius=0)])
	ellipsepoints_out = np.vstack([p for p in coords if not oval.contains_point(p, radius=0)])

	in_x = ellipsepoints_in[:,0]
	in_y = ellipsepoints_in[:,1]

	out_x = ellipsepoints_out[:,0]
	out_y = ellipsepoints_out[:,1]

	#----------------
	# NH3 Analysis


	print('\n(Over?)analysing NH3\nAssuming a cloud-forming altitude of {} mbar\n'.format(nh3_alt))

	idx_use = find_nearest(pressure, nh3_alt)


	hdul = fits.open(epochs[cccp] + '/' + nh3_file)
	nh3_vmr = hdul['VAL1'].data
	nh3_fsh = hdul['VAL2'].data
	hdul.close()

	nh3_vmr *= 1E+6 # ppm

	SVP = np.exp(A + np.divide(B, temp) + np.multiply(C, temp)) # bar
	SVP_show = SVP[idx_use] * 1e+6

	# The plot
	# Analysis

	for kk in range(2):

		if kk == 0:
			x = SVP_show[in_y, in_x].flatten()
			y = nh3_fsh[in_y, in_x].flatten()

		if kk == 1:
			x = SVP_show[out_y, out_x].flatten()
			y = nh3_fsh[out_y, out_x].flatten()


		cover = np.isfinite(x) & np.isfinite(y)

		x = x[cover]
		y = y[cover]

		c = np.polyfit(x, y, 1)
		p = np.poly1d(c)
		px = np.linspace(min(x), max(x), num=10)
		py = p(px)

		pearson = np.corrcoef(x, y)

		ax[cccp, 0].scatter(x, y, s=1, c=colours[kk])
		ax[cccp, 0].plot(px, py, color=colours[kk], linestyle='-', linewidth=1.0, label=sub_labels[kk] + ' p={aa:.2f}'.format(aa=pearson[0][1]))

	ax[cccp, 0].set_xlabel('{aa} mbar'.format(aa=nh3_alt) + ' SVP (x10$^{-6}$ bar)')
	ax[cccp, 0].set_ylabel('NH$_{3}$ FSH')
	ax[cccp, 0].set_title(r"$\bf{" + labels[cccp][0] + "}$" + " " + epochs[cccp].replace('j', 'J').replace('aug', 'August') + " NH$_3$ analysis", loc='left')

	ax[cccp, 0].grid()
	ax[cccp, 0].legend()


	#----------------
	# PH3 Analysis


	print('\n(Over?)analysing PH3\n')


	hdul = fits.open(epochs[cccp] + '/' + ph3_file)
	ph3_vmr = hdul['VAL1'].data
	ph3_fsh = hdul['VAL2'].data
	hdul.close()

	ph3_vmr *= 1E+6 # ppm

	hdul = fits.open(epochs[cccp] + '/' + aero_file)
	op = hdul['VAL2'].data
	hdul.close()

	# The second layer

	for kk in range(2):

		if kk == 0:
			x = op[in_y, in_x].flatten()
			y = ph3_fsh[in_y, in_x].flatten()

		if kk == 1:
			x = op[out_y, out_x].flatten()
			y = ph3_fsh[out_y, out_x].flatten()

		cover = np.isfinite(x) & np.isfinite(y)
		x = x[cover]
		y = y[cover]

		c = np.polyfit(x, y, 1)
		p = np.poly1d(c)
		px = np.linspace(min(x), max(x), num=10)
		py = p(px)

		pearson = np.corrcoef(x, y)

		ax[cccp, 1].scatter(x, y, s=1, c=colours[kk])
		ax[cccp, 1].plot(px, py, color=colours[kk], linestyle='-', linewidth=1.0, label=sub_labels[kk] + ' p={aa:.2f}'.format(aa=pearson[0][1]))

	ax[cccp, 1].set_xlabel('Opacity')
	ax[cccp, 1].set_ylabel('PH$_3$ FSH')
	ax[cccp, 1].set_title(r"$\bf{" + labels[cccp][1] + "}$" + " " + epochs[cccp].replace('j', 'J').replace('aug', 'August') + " PH$_3$ analysis", loc='left')

	ax[cccp, 1].legend()
	ax[cccp, 1].grid()


fig.subplots_adjust(wspace=wspace, hspace=hspace)
plt.savefig(fig_name, bbox_inches='tight')
plt.show()


#===========================================================================

print('End of script\n')
