#==================================================================================
# Description
"""
Loads the bkg data and generates a median bkg spectrum, applies to the d_combined
data input below.

Usage
		[Input parameters below]

		python -W ignore apply_bkg.py


"""
#==================================================================================
# Inputs


bkg = 'july/bkg/stage3/d1' # location of bkg files

band = 'ch3-medium' # band to generate bkg for

final_dir = 'july/east/d_combined_nofringe' # data to apply bkg to


#==================================================================================
# Imports


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#==================================================================================
# Main


if not os.path.exists(final_dir + '_bkg'):
	os.mkdir(final_dir + '_bkg')


hdul = fits.open(bkg + '/Level3_' + band + '_s3d.fits')
bkg_cube = hdul['SCI'].data
hdul.close()

bkg_wave = gen_wave_arr(bkg + '/Level3_' + band + '_s3d.fits')

arral = np.where(bkg_cube == 0.0)
bkg_cube[arral] = np.nan

bkg_spec = np.nanmedian(np.nanmedian(bkg_cube, axis=1), axis=1)


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))

ax.plot(bkg_wave, bkg_spec, 'k-', linewidth=0.7)

ax.set_xlabel('Wavelength (µm)')
ax.set_ylabel('Surf Brightness (MJysr$^{-1}$)')

ax.grid()
plt.show()


hdul = fits.open(final_dir + '/Level3_' + band + '_s3d_nav.fits')

spec_data = hdul['SCI'].data

spec_data = np.subtract(spec_data, bkg_spec[:, None, None])

hdul['SCI'].data = spec_data
hdul.writeto(final_dir + '_bkg/Level3_' + band + '_s3d_nav.fits', overwrite=True)


#==================================================================================

print('End of script\n')
