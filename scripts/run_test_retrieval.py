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


bands = ['ch1-long', 'ch2-short']

wmin = [7.30, 7.58085]
wmax = [7.58, 8.10]

ret_dir = 'stage2'
prior_dir = 'stage1' # results from the previous retrieval
				# Only required if ret_dir != stage1, usually prior_dir == stage1
temp_file = '' # only required if stage == 1
xsc_file = 'grs_single.xsc'
aero_file = '' # Only required if stage == 1
mre_skip = 2053 # Line temperature data begins in the mre file, only needed if stage == 2

delta_spec = 8.0   # Value to multiply the uncertainties by

tp_cons = 0.1 # Value to multiply dT by below turn_alt (only effective if stage == 2)
turn_alt = 100 # mbvar, altitude to return to normal dT (if stage == 2)

constrainQ = True

walltime = '18:00:00'
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
# Functions


# Produces the spx file for each lon/lat
def gen_wavefiles(bands, lon, lat, p):

	wave = [None]*len_b
	spec_err = [None]*len_b

	for kk in range(len_b):
		path = ''

		if os.path.exists('stage0/' + bands[kk] + '/spx/{a}_{b}.spx'.format(a=lon, b=lat)):
			path = 'stage0/' + bands[kk] + '/spx/{a}_{b}.spx'.format(a=lon, b=lat)

		data1 = np.loadtxt(path, skiprows=4, usecols=(0, 1))
		spec_err[kk] = np.loadtxt(path.replace('spx/', 'fmerror/').replace('spx', 'txt'), usecols=(1))

		wave[kk] = data1[:,0]
		spec = data1[:,1]

		idx_start = find_nearest(wave[kk], wmin[kk])
		idx_stop = find_nearest(wave[kk], wmax[kk])

		spec = spec[idx_start:idx_stop+1]
		spec_err[kk] = spec_err[kk][idx_start:idx_stop+1]
		wave[kk] = wave[kk][idx_start:idx_stop+1]

		# Special Wavegrid Instructions

		if bands[kk] == 'ch2-short' and constrainQ == True:
			print('Dividing Q branch uncertainties by 10\n')

			arral = np.where(np.logical_and(wave[kk] > 7.6560, wave[kk] < 7.67821))
			spec_err[kk][arral] /= 10.0

			for k in range(10):
				spec_err[kk][max(arral[0])+ (k+1)] = spec_err[kk][max(arral[0])+ (k+1)]/(10-k)
				spec_err[kk][arral[0][0] - (10-k)] = spec_err[kk][arral[0][0] - (10-k)]/(k+1)


		with open(path, 'r') as g:
			data = g.readlines()

		data_array = data[3]
		g.close()

		len_spec = len(spec)


		if not os.path.exists(ret_dir + '/core{}/jwst_mrs.spx'.format(p)):
			h = open(ret_dir + '/core{}/jwst_mrs.spx'.format(p), 'a')
			h.write('0.0000 {a} {b} {c}\n'.format(a=lat, b=lon, c=len_b))

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

	np.savetxt(ret_dir + '/core{}/fmerror_data.txt'.format(p), b, fmt='%.6e', delimiter='	')


#=======================================================
# Main


if not os.path.exists(ret_dir):
	os.mkdir(ret_dir)

len_b = len(bands)

if 'stage1' in ret_dir or 'stage2_aux' in ret_dir:
	print('\nInitilising stage1 or stage2_aux\n')

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
	print('\nInitilising plain stage 2\n')

	# Loading the ID parameters from the previous retrieval step
	old_id_data = np.loadtxt(prior_dir + '/id_conserve.txt', skiprows=1)
	old_id = old_id_data[:,0]
	lon_arr = old_id_data[:,1]
	lat_arr = old_id_data[:,2]
	len_old_id = len(old_id)

	pressure = np.loadtxt('jupiter_v2022miri.ref', skiprows=25, usecols=(1)) # atm


if os.path.exists(ret_dir + '/id_conserve.txt'):
	os.remove(ret_dir + '/id_conserve.txt')

f = open(ret_dir + '/id_conserve.txt', 'a')

f.write("#ID	lon	lat\n")

n = 1


if 'stage1' in ret_dir or 'stage2_aux' in ret_dir:
	for j in range(len(lon_arr)):
		for i in range(len(lat_arr)):
			print('({aa}, {bb})'.format(aa=lon_arr[j], bb=lat_arr[i]))

			try:

				if os.path.exists(ret_dir + '/core{aa}'.format(aa=n)):
					shutil.rmtree(ret_dir + '/core{aa}'.format(aa=n))
				os.mkdir(ret_dir + '/core{aa}'.format(aa=n))

				gen_wavefiles(bands, lon_arr[j], lat_arr[i], n)

				shutil.copy(temp_file, ret_dir + '/core{}/tempapr.dat'.format(n))
				shutil.copy(aero_file, ret_dir + '/core{}/aerosol.ref'.format(n))
				shutil.copy('jupiter_v2022miri.ref', ret_dir + '/core{}/jupiter_v2022miri.ref'.format(n))
				shutil.copy('jupiter_v2022miri.prf', ret_dir + '/core{}/jupiter_v2022miri.prf'.format(n))
				shutil.copy(xsc_file, ret_dir + '/core{}/grs_single.xsc'.format(n))

				f.write("{aa}   {bb}    {cc}\n".format(aa=n, bb=lon_arr[i], cc=lat_arr[i]))
				n = n + 1

			except:
				print('No data found')
				shutil.rmtree(ret_dir + '/core{aa}'.format(aa=n))



if ret_dir == 'stage2':
	for i in range(len_old_id):

		print('({aa}, {bb})'.format(aa=lon_arr[i], bb=lat_arr[i]))

		try:

			if os.path.exists(ret_dir + '/core{aa}'.format(aa=n)):
				shutil.rmtree(ret_dir + '/core{aa}'.format(aa=n))
			os.mkdir(ret_dir + '/core{aa}'.format(aa=n))

			gen_wavefiles(bands, lon_arr[i], lat_arr[i], n)


			# Special Temperature Instructions

			temp_data = np.loadtxt(prior_dir + '/core{}/core_1/nemesis.mre'.format(int(old_id[i])), skiprows=mre_skip, max_rows=120, usecols=(4, 5))
			temp = temp_data[:,0]
			dtemp = temp_data[:,1]


			arral = np.where(pressure > turn_alt/1000)

			dtemp[arral] = dtemp[arral] * tp_cons

			idx_start = arral[0][-1]
			idx_finish = arral[0][-1] + 10

			n_points = 10
			dt = (dtemp[idx_finish] - dtemp[idx_start]) / n_points

			for k in range(10):
				dtemp[idx_start + (k + 1)] = dtemp[idx_start] + (dt*k)

			# Lower corner
			dtemp[idx_start+2] = np.median([dtemp[idx_start+1], dtemp[idx_start+3]])

			dtemp[idx_start+1] = np.median([dtemp[idx_start], dtemp[idx_start+2]])
			dtemp[idx_start+3] = np.median([dtemp[idx_start+2], dtemp[idx_start+4]])

			dtemp[idx_start] = np.median([dtemp[idx_start-1], dtemp[idx_start+1]])
			dtemp[idx_start+4] = np.median([dtemp[idx_start+3], dtemp[idx_start+5]])

			dtemp[idx_start-1] = np.median([dtemp[idx_start-2], dtemp[idx_start]])
			dtemp[idx_start+5] = np.median([dtemp[idx_start+4], dtemp[idx_start+6]])

			dtemp[idx_start-2] = np.median([dtemp[idx_start-3], dtemp[idx_start-1]])
			dtemp[idx_start+6] = np.median([dtemp[idx_start+5], dtemp[idx_start+7]])

			dtemp[idx_start-3] = np.median([dtemp[idx_start-4], dtemp[idx_start-2]])
			dtemp[idx_start+7] = np.median([dtemp[idx_start+6], dtemp[idx_start+8]])

			dtemp[idx_start+2] = dtemp[idx_start+2] + (((dtemp[idx_finish] - dtemp[idx_start])/6.0) * 0.02)

			# Upper corner
			dtemp[idx_finish] = np.median([dtemp[idx_finish-1], dtemp[idx_finish+1]])

			dtemp[idx_finish-1] = np.median([dtemp[idx_finish-2], dtemp[idx_finish]])
			dtemp[idx_finish+1] = np.median([dtemp[idx_finish], dtemp[idx_finish+2]])

			dtemp[idx_finish-2] = np.median([dtemp[idx_finish-3], dtemp[idx_finish-1]])
			dtemp[idx_finish+2] = np.median([dtemp[idx_finish+1], dtemp[idx_finish+3]])

			dtemp[idx_finish-3] = np.median([dtemp[idx_finish-4], dtemp[idx_finish-2]])
			dtemp[idx_finish+3] = np.median([dtemp[idx_finish+2], dtemp[idx_finish+4]])

			dtemp[idx_finish-4] = np.median([dtemp[idx_finish-5], dtemp[idx_finish-3]])
			dtemp[idx_finish+4] = np.median([dtemp[idx_finish+3], dtemp[idx_finish+5]])

			dtemp[idx_finish-5] = np.median([dtemp[idx_finish-6], dtemp[idx_finish-4]])
			dtemp[idx_finish+5] = np.median([dtemp[idx_finish+4], dtemp[idx_finish+6]])

			dtemp[idx_finish] = dtemp[idx_finish] - (((dtemp[idx_finish] - dtemp[idx_start])/6.0) * 0.02)


			a = [None]*3
			a[0] = pressure
			a[1] = temp
			a[2] = dtemp

			b = np.transpose(a)
			np.savetxt(ret_dir + '/core{}/tempapr.dat'.format(n), b, fmt='%.5e', delimiter='	', header='120	1.50000', comments='')

			shutil.copy(prior_dir + '/core{}/core_1/aerosol.prf'.format(int(old_id[i])), ret_dir + '/core{}/aerosol.ref'.format(n))
			shutil.copy(prior_dir + '/core{}/core_1/nemesis.prf'.format(int(old_id[i])), ret_dir + '/core{}/prior.ref'.format(n))
			shutil.copy(prior_dir + '/core{}/core_1/nemesis.prf'.format(int(old_id[i])), ret_dir + '/core{}/prior.prf'.format(n))
			shutil.copy(xsc_file, ret_dir + '/core{}/grs_single.xsc'.format(n))
			shutil.copy('preNemesis.py', ret_dir + '/core{}/preNemesis.py'.format(n))

			with open(ret_dir + '/core{}/prior.ref'.format(n), 'r') as fg:
				datafg = fg.readlines()
			fg.close()

			datafg[0] = "           1\n           1\n"

			with open(ret_dir + '/core{}/prior.ref'.format(n), 'w') as fg:
				fg.writelines(datafg)

			os.chdir(ret_dir + '/core{}'.format(n))
			call('python preNemesis.py',shell=True)
			os.chdir('../..')

			f.write("{aa}	{bb}	{cc}\n".format(aa=n, bb=lon_arr[i], cc=lat_arr[i]))
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
