#=======================================================================
"""
Extracts spectra from a specified cube, converts to spx file units
(W/cm2/sr/um) and saves in a spx file.
Also a corresponding fmerror file is generated for each spectrum

Usage:
		[Input parameters below]

		python retrieval_map_multibands.py


"""
#=======================================================================
# Inputs


working_dir = 'july/centre' # Structure should be: [epoch]/[tile]

stage = 0
band_name = 'ch1-long' # Only needed if stage == 0

dir_in = 'd_combined_0.5' # Location of Nyquist dither-combined data

mult_err = 1.0 # Value to multiply spectral uncertainty by (see note below)

remove_pix = False # Set to True to remove spaxels
positions = []


#--------------------
"""
Note about stages

As of 04/04/23:

Stage 1 is the stratospheric temperature retrieval, every wavelength in
the range 7.30 µm - 8.10 µm is used. The spectral uncertainties are
multiplied by 50, the uncertainties around the Q-band are only multiplied
by 5

Stage 2 is the everything retrieval, it uses the results of stage 1 as the
prior to retrieve everything else. Every other spectral point is taken in the
range 6.0 µm (start of ch1-long) - 10.75 µm. The spectral uncertainties are
multiplied by 50 except for the Q-branch which are multiplied by 5.

Stage 23 is the alternate spectral range (8.1 - 10.75 µm) that can be used in
stage 2, in this case every spectral point will be included

Enter stage above
"""
#--------------------


#======================================================================
# Imports


import numpy as np
from astropy.io import fits
import os
import shutil
import glob

c = 2.99792458e+8
h=6.626e-34
kb=1.3806e-23
c1=2*h*c*c
c2=(h*c)/kb


#=====================================================================
# Functions


#Finds nearest array wavelength value to specified range
#Returns index of that wavelength
def find_nearest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx


# Converts Jy/sr to spx file units of W/cm2/sr/µm
def convert_surfb_spx(rad, wave):
	# rad starts in units of Jy/sr

	l = wave * (1e-6)

	rad = rad * (1e-26) # W/m2/sr/Hz
	rad = rad * (c / l**2) # W/m2/sr/m
	rad = rad / (1e4) # W/cm2/sr/m
	rad = rad / (1e6) # W/cm2/sr/µm

	return rad


#=======================================================================
#Main


# Setting up the retrieval parameters depending on stage selected
if stage == 0:
	bands = [band_name]


if stage == 1 or stage == 11:
	bands = ['ch1-long', 'ch2-short', 'ch2-medium', 'ch2-long']

	wmin = [7.3007, 7.58, 8.72, 10.08]
	wmax = [7.58, 8.72, 10.08, 10.75]

if stage == 2 or stage == 24 or stage == 23 or stage == 22:
	bands = ['ch2-short', 'ch2-medium', 'ch2-long']

	wmin = [8.10, 8.72, 10.08]
	wmax = [8.72, 10.08, 10.75]

if stage == 3:
	if working_dir == 'aug/east':
		print('\nProcessing aug/east\nAlternate method being used\n')

		bands = ['ch2-short']
		wmin = [7.50]
		wmax = [8.80]

	else:
		print('\nStandard method being used\n')

		bands = ['ch1-long', 'ch2-short']

		wmin = [7.30, 7.58]
		wmax = [7.58, 8.10]

len_b = len(bands)

if not os.path.exists(working_dir + '/retrieval_data'):
	os.mkdir(working_dir + '/retrieval_data')

if not os.path.exists(working_dir + '/retrieval_data/stage{}'.format(stage)):
	os.mkdir(working_dir + '/retrieval_data/stage{}'.format(stage))

if stage == 0:
	if os.path.exists(working_dir + '/retrieval_data/stage{}/'.format(stage) + band_name):
		shutil.rmtree(working_dir + '/retrieval_data/stage{}/'.format(stage) + band_name)

	os.mkdir(working_dir + '/retrieval_data/stage{}/'.format(stage) + band_name)
	os.mkdir(working_dir + '/retrieval_data/stage{}/'.format(stage) + band_name + '/spx')
	os.mkdir(working_dir + '/retrieval_data/stage{}/'.format(stage) + band_name + '/fmerror')

if not stage == 0:
	os.mkdir(working_dir + '/retrieval_data/stage{}/spx'.format(stage))
	os.mkdir(working_dir + '/retrieval_data/stage{}/fmerror'.format(stage))

spec_data = [None]*len_b
spec_err = [None]*len_b

sa = [None]*len_b
ea = [None]*len_b
aa = [None]*len_b

lon = [None]*len_b
lat = [None]*len_b

wave = [None]*len_b

dx = [None]*len_b
dy = [None]*len_b
dz = [None]*len_b

#Finds and opens cubes
for kk in range(len_b):
	print(bands[kk])

	hdul = fits.open(working_dir + '/' + dir_in + '/Level3_' + bands[kk] + '_s3d_nav.fits')

	spec_data[kk] = hdul['SCI'].data  # MJy/sr
	spec_err[kk] = hdul['ERR'].data # MJy/sr

	spec_err[kk] *= mult_err

	sa[kk] = hdul['SOLAR_ANGLE'].data
	ea[kk] = hdul['EMISSION_ANGLE'].data
	aa[kk] = hdul['AZIMUTH_ANGLE'].data

	lon[kk] = hdul['LON_WEST'].data
	lat[kk] = hdul['LAT_PGR'].data

	scihdr = hdul['SCI'].header

	wave[kk] = np.arange(scihdr['NAXIS3'])*scihdr['CDELT3']+scihdr['CRVAL3']

	# Special instructions

	if remove_pix == True:
		print('\nRemoving Positions\n')

		spec_data[kk] = np.delete(spec_data[kk], positions, axis=0)
		spec_err[kk] = np.delete(spec_err[kk], positions, axis=0)
		wave[kk] = np.delete(wave[kk], positions, axis=0)

	hdul.close()

	if not stage == 0:
		idx_start = find_nearest(wave[kk], wmin[kk])
		idx_stop = find_nearest(wave[kk], wmax[kk])

		if wave[kk][idx_start] < wmin[kk]:
			idx_start = idx_start + 1
		if wave[kk][idx_stop] > wmax[kk]:
			idx_stop = idx_stop - 1

		spec_data[kk] = spec_data[kk][idx_start:idx_stop+1]
		spec_err[kk] = spec_err[kk][idx_start:idx_stop+1]

		wave[kk] = wave[kk][idx_start:idx_stop+1]

		if stage == 1 or stage == 24:
			wave[kk] = wave[kk][::4]
			spec_data[kk] = spec_data[kk][::4]
			spec_err[kk] = spec_err[kk][::4]

		if stage == 11:
			if kk == 0 or kk == 2 or kk == 3:
				wave[kk] = wave[kk][::4]
				spec_data[kk] = spec_data[kk][::4]
				spec_err[kk] = spec_err[kk][::4]

			if kk == 1:
				idx1 = find_nearest(wave[kk], 8.1)
				idx2 = find_nearest(wave[kk], 8.7)

				wave1 = wave[kk][0:idx1][::4]
				wave2 = wave[kk][idx1:idx2]
				wave3 = wave[kk][idx2:-1][::4]
				wave[kk] = np.concatenate([wave1, wave2, wave3])

				spec_data1 = spec_data[kk][0:idx1][::4]
				spec_data2 = spec_data[kk][idx1:idx2]
				spec_data3 = spec_data[kk][idx2:-1][::4]
				spec_data[kk] = np.concatenate([spec_data1, spec_data2, spec_data3], axis=0)

				spec_err1 = spec_err[kk][0:idx1][::4]
				spec_err2 = spec_err[kk][idx1:idx2]
				spec_err3 = spec_err[kk][idx2:-1][::4]
				spec_err[kk] = np.concatenate([spec_err1, spec_err2, spec_err3], axis=0)


		if stage == 23:
			wave[kk] = wave[kk][::3]
			spec_data[kk] = spec_data[kk][::3]
			spec_err[kk] = spec_err[kk][::3]

		if stage == 22:
			wave[kk] = wave[kk][::2]
			spec_data[kk] = spec_data[kk][::2]
			spec_err[kk] = spec_err[kk][::2]



	dz[kk], dy[kk], dx[kk] = spec_data[kk].shape

	print(dx[kk], dy[kk], dz[kk])

	print('\n')



spec_data_noblank = [None]*len_b
spec_data_noblankt = [None]*len_b

spec_err_noblank = [None]*len_b
spec_err_noblankt = [None]*len_b

sa_noblankt = [None]*len_b
ea_noblankt = [None]*len_b
aa_noblankt = [None]*len_b

lon_noblankt = [None]*len_b
lat_noblankt = [None]*len_b

print(np.array(spec_data).shape)

if working_dir == 'aug/east' or stage == 0:
	print('Using aug/east or stage 0 arral approach\n')

	arral = np.where(np.isnan(spec_data[0][0]) == False)

else:
	print('Using other tile arral approach\n')

	hdul2 = fits.open(working_dir + '/' + dir_in + '/Level3_ch1-long_s3d_nav.fits')
	spec_data2 = hdul2['SCI'].data
	scihdr2 = hdul2['SCI'].header
	hdul2.close()

	hdul3 = fits.open(working_dir + '/' + dir_in + '/Level3_ch2-short_s3d_nav.fits')
	spec_data3 = hdul3['SCI'].data
	scihdr3 = hdul3['SCI'].header
	hdul3.close()

	wave2 = np.arange(scihdr2['NAXIS3'])*scihdr2['CDELT3']+scihdr2['CRVAL3']
	wave3 = np.arange(scihdr3['NAXIS3'])*scihdr3['CDELT3']+scihdr3['CRVAL3']

	arral = np.where(np.logical_and(np.isnan(spec_data2[find_nearest(wave2, 7.30)]) == False, np.isnan(spec_data3[find_nearest(wave3, 7.58)]) == False))

for kk in range(len_b):

	spec_data_noblank[kk] = [None]*dz[kk]
	spec_err_noblank[kk] = [None]*dz[kk]

	for k in range(dz[kk]):
		spec_data_noblank[kk][k] = spec_data[kk][k][arral]
		spec_err_noblank[kk][k] = spec_err[kk][k][arral]

	spec_data_noblankt[kk] = np.transpose(spec_data_noblank[kk])
	spec_err_noblankt[kk] = np.transpose(spec_err_noblank[kk])

	sa_noblankt[kk] = sa[kk][arral]
	ea_noblankt[kk] = ea[kk][arral]
	aa_noblankt[kk] = aa[kk][arral]

	lon_noblankt[kk] = lon[kk][arral]
	lat_noblankt[kk] = lat[kk][arral]

id = np.arange(1,(dy[0]*dx[0])+1)

d = np.reshape(id, (dy[0], dx[0]))
hdu2 = fits.PrimaryHDU(d)
hdu2.writeto(working_dir + '/retrieval_data/stage{}/retrieval_map.fits'.format(stage), overwrite=True)

id2 = d[arral]
len_id2 = len(id2)
id2_number = np.arange(1, len_id2+1)

f = [None]*2
f[0] = id2_number
f[1] = id2

g = np.transpose(f)

np.savetxt(working_dir + '/retrieval_data/stage{}/id_conserve.txt'.format(stage), g, fmt='%.2f', delimiter='	', header='Pixel	Position')

flux = [None]*len_b
flux_err = [None]*len_b

for kk in range(len_b):

	spec_data_noblankt[kk] *= (1e+6) # Jy/sr
	spec_err_noblankt[kk] *= (1e+6) # Jy/sr

	flux[kk] = convert_surfb_spx(spec_data_noblankt[kk], wave[kk])
	flux_err[kk] = convert_surfb_spx(spec_err_noblankt[kk], wave[kk])


for jj in range(len_id2):
	if not stage == 0:
		f = open(working_dir + '/retrieval_data/stage{a}/spx/{b}.spx'.format(a=stage, b=jj+1), 'a')
	if stage == 0:
		f = open(working_dir + '/retrieval_data/stage0/' + band_name + '/spx/{a}_{b}.spx'.format(a=lon_noblankt[0][jj], b=lat_noblankt[0][jj]), 'a')

	f.write('0.0000 {a} {b} {c}\n'.format(a=lat_noblankt[0][jj], b=lon_noblankt[0][jj], c=len_b))

	fm_wave = []
	fm_err = []

	for kk in range(len_b):

		f.write('	{d}\n       1\n'.format(d=len(wave[kk])))

		f.write('{a}	  {b}	   {e:.4f}	{f:.4f}      {g:.4f}  1.0000\n'.format(a=lat_noblankt[0][jj], b=lon_noblankt[0][jj], e=sa_noblankt[kk][jj], f=ea_noblankt[kk][jj], g=aa_noblankt[kk][jj]))

		an = [None]*3
		an[0] = wave[kk]
		an[1] = flux[kk][jj]
		#an[2] = [0.000000e+00]*dz[kk]
		an[2] = flux_err[kk][jj] / mult_err

		bn = np.transpose(an)

		np.savetxt(f, bn, fmt='%.6e', delimiter='	')

		f.write('\n')

		fm_wave.append(wave[kk])
		fm_err.append(flux_err[kk][jj])

	fm_wave = np.concatenate((fm_wave))
	fm_err = np.concatenate((fm_err))

	if stage == 3:
		arral = np.where(np.logical_and(fm_wave > 7.6560, fm_wave < 7.67821))
		fm_err[arral] /= 10.0

		for k in range(10):
			fm_err[max(arral[0])+ (k+1)] = fm_err[max(arral[0])+ (k+1)]/(10-k)
			fm_err[arral[0][0] - (10-k)] = fm_err[arral[0][0] - (10-k)]/(k+1)

	f.close()

	ann = [None]*2
	ann[0] = fm_wave
	ann[1] = fm_err
	bnn = np.transpose(ann)

	if not stage == 0:
		np.savetxt(working_dir + '/retrieval_data/stage{a}/fmerror/{b}.txt'.format(a=stage, b=jj+1), bnn, fmt='%.6e', delimiter='	')
	if stage == 0:
		np.savetxt(working_dir + '/retrieval_data/stage0/' + band_name + '/fmerror/{a}_{b}.txt'.format(a=lon_noblankt[0][jj], b=lat_noblankt[0][jj]), bnn, fmt='%.6e', delimiter=' 	')


#=======================================================================


print('End of script\n')
