#==============================================================================
# Description
"""
Compares the spectra and uncertainties from files input below
Use to determine what factor to multiply the new fmerror files by in the
all-up testing process

Usage

		[Input parameters below]

		python compare_uncertainties.py


"""
#==============================================================================
# Inputs


old_file = '/data/nemesis/jwst/MIRI_IFU/new_pipeline/JupiterGRS_2022/self_flat/july/centre/\
stage5_despike/d1_fringe_nav/Level3_ch2-medium_s3d_nav.fits'

new_file = 'july/centre/stage5_despike/d1_fringe/Level3_ch2-medium_s3d_nav.fits'

x_pos = 10 # 1 idx
y_pos = 14 # 1 idx

fig_width = 12
fig_height = 7


#==============================================================================
# Imports


import numpy as np
import sys
from astropy.io import fits

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#==============================================================================
# Main


#fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(fig_width, fig_height))

colours = ['black', 'red']
labels = ['Old file', 'New file']

for kk in range(2):
	print('\n')

	if kk == 0:
		file = old_file
	if kk == 1:
		file = new_file

	print(labels[kk])

	hdul = fits.open(file)
	spec_grid = hdul['SCI'].data
	dspec_grid = hdul['ERR'].data
	lon_grid = hdul['LON'].data
	lat_grid = hdul['LAT_GRAPHIC'].data

	scihdr = hdul['SCI'].header
	wave = np.arange(scihdr['NAXIS3'])*scihdr['CDELT3'] + scihdr['CRVAL3']
	hdul.close()

	print('Lon/Lat position = ({aa}, {bb})'.format(aa=lon_grid[y_pos-1, x_pos-1], bb=lat_grid[y_pos-1, x_pos-1]))

	spec = spec_grid[:, y_pos-1, x_pos-1]

	if kk == 0:
		dspec_old = dspec_grid
	if kk == 1:
		dspec_new = dspec_grid


print(np.nanmedian(dspec_new/dspec_old))

"""
	ax[0].plot(wave, spec, linestyle='-', linewidth=0.5, color=colours[kk], label=labels[kk])
	#ax[0].fill_between(wave, spec-dspec, spec+dspec, color=colours[kk], alpha=0.3)


# Plot 1
ax[0].set_xlabel('Wavelength (µm)')
ax[0].set_ylabel('Brightness (MJysr$^{-1}$)')
ax[0].set_title('Spectra')

ax[0].grid()
ax[0].legend()


# Plot 2

print('Median difference = {}'.format(np.nanmedian(dspec_new/dspec_old)))

ax[1].plot(wave, dspec_new/dspec_old, linewidth=0.5, color='black')
ax[1].set_xlabel('Wavelength (µm)')
ax[1].set_ylabel('dBrightness (MJysr$^{-1}$)')
ax[1].set_title('Factorial difference in uncertainties')

ax[1].grid()


# Final config

plt.suptitle('Pix position = ({aa}, {bb})'.format(aa=x_pos, bb=y_pos))
plt.show()
"""


#==============================================================================

print('\nEnd of script\n')
