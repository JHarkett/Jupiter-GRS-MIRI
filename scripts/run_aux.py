#=======================================================
"""
Detects spx files input below and organises into 
directories ready to be run using NEMESIS


Usage

		[Optional: Edit special instructions starting on line 179]

		[Input parameters below]

		python run_test_retrieval.py


"""
#=======================================================
# Inputs


bands = ['ch2-short', 'ch2-medium', 'ch2-long']

wmin = [8.10, 8.72, 10.08]
wmax = [8.72, 10.08, 10.75]

ret_dir = 'stage1'
prior_dir = '' # results from the previous retrieval
				# Only required if ret_dir != stage1

temp_file = 'tempapr.dat' # only required if ret_dir == stage1
xsc_file = 'grs_single.xsc'
aero_file = 'aerosol.ref'

delta_spec = 5.0   # Value to multiply the uncertainties by

walltime = '22:00:00'
vmem = '16G'


#=======================================================
# Imports


import numpy as np
import os
import shutil
from subprocess import call
import glob
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=======================================================
# Main


if not os.path.exists(ret_dir):
	os.mkdir(ret_dir)

len_b = len(bands)

if 'stage1' in ret_dir:
	print('\nInitilising stage1\n')

	min_lon = []
	max_lon = []
	min_lat = []
	max_lat = []

	for jj in range(len_b):

		sstring = 'stage0/' + bands[jj] + '/spx/*.spx'

		spx_files = sorted(glob.glob(sstring))

		if len(spx_files) != 0:

			min_lon.append(float(spx_files[0].replace('stage0/' + bands[jj] + '/spx/', '').replace('.spx', '').replace('_', ' ').split()[0]))
			max_lon.append(float(spx_files[-1].replace('stage0/' + bands[jj] + '/spx/', '').replace('.spx', '').replace('_', ' ').split()[0]))

			sstring = 'stage0/' + bands[jj] + '/spx/{}_*.spx'.format(min_lon[-1]+5.0)
			spx_files = sorted(glob.glob(sstring))

			max_lat.append(float(spx_files[0].replace('stage0/' + bands[jj] + '/spx/{}_'.format(min_lon[-1]+5.0), '').replace('.spx', '')))
			min_lat.append(float(spx_files[-1].replace('stage0/' + bands[jj] + '/spx/{}_'.format(min_lon[-1]+5.0), '').replace('.spx', '')))

	min_lon_total = np.amin(min_lon) - 1.0
	max_lon_total = np.amax(max_lon) + 1.0

	min_lat_total = np.amin(min_lat) - 1.0
	max_lat_total = np.amax(max_lat) + 1.0

	print('\nLongitude range: {aa} - {bb} degrees'.format(aa=max_lon_total, bb=min_lon_total))
	print('Latitude range: {aa} - {bb} degrees\n'.format(aa=min_lat_total, bb=max_lat_total))


	lon_arr = np.arange(max_lon_total, min_lon_total-0.1, -0.5)
	lat_arr = np.arange(min_lat_total, max_lat_total+0.1, 0.5)

else:
	print('\nInitilising stage 2\n')

	hdul = fits.open(prior_dir + '/lon_lat.fits')
	lon_grid = hdul['LON_WEST'].data
	lat_grid = hdul['LAT_PGR'].data
	hdul.close()

	lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

	pressure = np.loadtxt('jupiter_v2022miri.ref', skiprows=25, usecols=(1)) # atm

	hdul = fits.open(prior_dir + '/000_retrieved_cleaned.fits')
	temp = hdul['VARIABLE'].data
	temp_err = hdul['ERR'].data
	hdul.close()


if os.path.exists(ret_dir + '/id_conserve.txt'):
	os.remove(ret_dir + '/id_conserve.txt')

f = open(ret_dir + '/id_conserve.txt', 'a')

f.write("#ID	lon	lat\n")

n = 1


for j in range(len(lon_arr)):
	for i in range(len(lat_arr)):
		print('({aa}, {bb})'.format(aa=lon_arr[j], bb=lat_arr[i]))

		if n > 4:
			raise Exception('STOP!')

		try:

			if os.path.exists(ret_dir + '/core{aa}'.format(aa=n)):
				shutil.rmtree(ret_dir + '/core{aa}'.format(aa=n))
			os.mkdir(ret_dir + '/core{aa}'.format(aa=n))

			wave = [None]*len_b
			spec_err = [None]*len_b

			for kk in range(len_b):
				path = ''

				if os.path.exists('stage0/' + bands[kk] + '/spx/{a}_{b}.spx'.format(a=lon_arr[j], b=lat_arr[i])):
					path = 'stage0/' + bands[kk] + '/spx/{a}_{b}.spx'.format(a=lon_arr[j], b=lat_arr[i])

				data1 = np.loadtxt(path, skiprows=4, usecols=(0, 1))
				spec_err[kk] = np.loadtxt(path.replace('spx/', 'fmerror/').replace('spx', 'txt'), usecols=(1))

				wave[kk] = data1[:,0]
				spec = data1[:,1]

				idx_start = find_nearest(wave[kk], wmin[kk])
				idx_stop = find_nearest(wave[kk], wmax[kk])

				spec = spec[idx_start:idx_stop+1]
				spec_err[kk] = spec_err[kk][idx_start:idx_stop+1]

				wave[kk] = wave[kk][idx_start:idx_stop+1]

				with open(path, 'r') as g:
					data = g.readlines()

				data_array = data[3]
				g.close()

				len_spec = len(spec)


				if not os.path.exists(ret_dir + '/core{}/jwst_mrs.spx'.format(n)):
					h = open(ret_dir + '/core{}/jwst_mrs.spx'.format(n), 'a')
					h.write('0.0000 {a} {b} {c}\n'.format(a=lat_arr[i], b=lon_arr[j], c=len_b))

				h.write('      {}\n	  1\n'.format(len_spec))
				h.write(data_array)

				for k in range(len_spec):
					h.write('{a:.6e}	{b:.6e}	{c:.6e}\n'.format(a=wave[kk][k], b=spec[k], c=spec_err[kk][k]))

			h.close()

			wave = np.concatenate(wave)
			spec_err = np.concatenate(spec_err) * delta_spec

			a = [None]*2
			a[0] = wave
			a[1] = spec_err

			b = np.transpose(a)

			np.savetxt(ret_dir + '/core{}/fmerror_data.txt'.format(n), b, fmt='%.6e', delimiter='	')

			if 'stage1' in ret_dir:
				shutil.copy(temp_file, ret_dir + '/core{}/tempapr.dat'.format(n))

			else:
				a = [None]*3
				a[0] = pressure
				a[1] = temp[:, i, j]
				a[2] = temp_err[:, i, j]

				b = np.transpose(a)
				np.savetxt(ret_dir + '/core{}/tempapr.dat'.format(n), b, fmt='%.5e', delimiter='	', header='120	1.50000', comments='')

			shutil.copy('preNemesis.py', ret_dir + '/core{}/preNemesis.py'.format(n))
			shutil.copy(aero_file, ret_dir + '/core{}/aerosol.ref'.format(n))
			shutil.copy('jupiter_v2022miri.ref', ret_dir + '/core{}/jupiter_v2022miri.ref'.format(n))
			shutil.copy('jupiter_v2022miri.prf', ret_dir + '/core{}/jupiter_v2022miri.prf'.format(n))
			shutil.copy(xsc_file, ret_dir + '/core{}/grs_single.xsc'.format(n))

			os.chdir(ret_dir + '/core{}'.format(n))
			call('python preNemesis.py',shell=True)
			os.chdir('../..')

			f.write("{aa}	{bb}	{cc}\n".format(aa=n, bb=lon_arr[j], cc=lat_arr[i]))
			n = n + 1

		except:
			print('No data found')
			shutil.rmtree(ret_dir + '/core{aa}'.format(aa=n))


f.close()


print('Generating qsub file\nwalltime = {pp}\nvmem = {pn}\n'.format(pp=walltime, pn=vmem))

ids = np.loadtxt(ret_dir + '/id_conserve.txt', skiprows=1, usecols=(0))
len_i = len(ids)

if os.path.exists(ret_dir + '/submitjob'):
        os.remove(ret_dir + '/submitjob')

f = open(ret_dir + '/submitjob','a')
f.write('#!/bin/bash\n')
f.write('#\n')
f.write('#SBATCH --job-name=nemesis\n')
f.write('#SBATCH --account=nemesis\n')
f.write('#SBATCH --time=' + walltime + '\n')
f.write('#SBATCH --mem=' + vmem + '\n')
f.write('#SBATCH --export=NONE\n')
f.write('#SBATCH --cpus-per-task=1\n')
f.write('#SBATCH --array=1-{}\n'.format(len_i))
f.write('export PATH=~/bin/ifort:$PATH\n')
f.write('cd $SLURM_SUBMIT_DIR\n')
f.write('inputdir=core${SLURM_ARRAY_TASK_ID}/core_1\n')
f.write('cd $inputdir\n')
f.write('Nemesis < nemesis.nam > log_${SLURM_ARRAY_TASK_ID}\n')
f.close()


#=======================================================

print('End of script\n')
