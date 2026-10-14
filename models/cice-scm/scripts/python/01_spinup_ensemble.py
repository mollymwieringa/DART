###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import sys
import xarray as xr
import numpy as np
import pandas as pd 

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the case name 
case_name = sys.argv[1]
user = sys.argv[2]
location = sys.argv[3]

# set the relevant directory names 
# where is this project repo located on your system?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# where does your forcing data live?
forcing_dir = project_dir + 'data/forcings/'
# where is the Icepack installation located on your system?
icepack_dir = '/glade/work/'+user+'/Icepack/'
# where is the scratch directory Icepack will write output to?
scratch_dir = '/glade/derecho/scratch/'+user+'/'
# where are you storing output
storage_dir = scratch_dir + '/ICEPACK_RUNS/' + case_name
# where is your template Icepack namelist file
namelist_template = icepack_dir + case_name + '/icepack_in'

# set machine name and compiler
machine = 'derecho'
compiler = 'intel'

# choose a spinup length and ensemble size (limited at this 
# time to 30 or less)
spinup_length = 11
ensemble_size = 30

if location != 'default':
    locs = {'Barents': [1.309, 0.698132], 
            'CoastalCanada': [1.41372, 6.24828], 
            'SibChuk': [1.3183906501, 3.0446500856], 
            'CentralArctic': [1.53589, 0]}

###############################################################
## BEGIN LAUNCH- DO NOT EDIT BELOW THIS LINE                 ##
###############################################################
#-----------------------------------------#
# 1. Set parameter details
#-----------------------------------------#
R_snw = 1.5
dt_mlt = 1.0

#-----------------------------------------#
# 2. Set up an Icepack case
#-----------------------------------------#
# go to the Icepack directory
os.chdir(icepack_dir)

# set up a new case 
comd = './icepack.setup -c '+case_name+' -m '+machine+' -e '+compiler
os.system(comd)

#-----------------------------------------#
# 3. Build the case 
#-----------------------------------------#
# go to the case directory
os.chdir(case_name)

# build the case 
comd = './icepack.build'
os.system(comd)

# check that the case was built correctly
if os.path.exists(storage_dir) is False:
    AssertionError('Model did not build correctly! Please rebuild model or check your directory setups in the Icepack installation.')
else:
    print('Model has been built!')

#-----------------------------------------#
# 4. Run the spinup
#-----------------------------------------#
spinup_year = 1

mem = 1
while mem <= ensemble_size:
    inst_string ='{0}'.format('%04d' % mem) 
    print('Running member '+inst_string+'...')

    # make tidy directories for model output
    os.chdir(storage_dir)
    os.makedirs('mem' + inst_string + '/history/')
    os.makedirs('mem' + inst_string + '/restart/')

    # link the model executable for each member to the main one built for the case
    os.symlink(storage_dir+'/icepack','mem' + inst_string+'/icepack')

    # begin working on an individual ensemble member
    os.chdir(storage_dir+'/mem' + inst_string)
    namelist = f90nml.read(namelist_template)

    # set case settings 
    namelist['setup_nml']['year_init'] = 2000
    namelist['setup_nml']['npt'] = 8760 * spinup_length
    namelist['setup_nml']['restart'] = False
    namelist['setup_nml']['runtype_startup'] = True
    namelist['setup_nml']['ice_ic'] = 'default'
    namelist['setup_nml']['dumpfreq'] = 'y'
    namelist['setup_nml']['restart_format'] = 'nc'
    namelist['setup_nml']['history_format'] = 'nc'
    namelist['setup_nml']['conserv_check'] = True
    if location != 'default':
        namelist['setup_nml']['input_lat'] = locs[location][0]
        namelist['setup_nml']['input_lon'] = locs[location][1]

    # change namelist parameters of note
    namelist['shortwave_nml']['r_snw'] = R_snw
    namelist['shortwave_nml']['dt_mlt'] = dt_mlt

    # set namelist forcing options
    namelist['forcing_nml']['data_dir'] = forcing_dir +location + '/spinup/'
    namelist['forcing_nml']['atm_data_type'] = 'JRA55'
    namelist['forcing_nml']['ocn_data_type'] = 'ISPOL'
    namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_'+inst_string+'.txt'
    namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_'+inst_string+'.txt'
    namelist['forcing_nml']['default_season'] = 'winter'
    namelist['forcing_nml']['lateral_flux_type'] = 'uniform_ice'
    namelist['forcing_nml']['update_ocn_f'] = True
    namelist['forcing_nml']['fyear_init'] = 2000
    namelist['forcing_nml']['ycycle'] = spinup_length
    # write namelist to needed file
    namelist.write('icepack_in',force=True)

    # run the member instance and dump output 
    comd = './icepack > icepack.out'
    os.system(comd)

    # advance to the next ensemble member 
    mem += 1

check_spinup = glob.glob(storage_dir + '/mem*/restart/iced.2011-01-01-00000.nc')

if len(check_spinup) != ensemble_size:
    print('HEY!! Spinup did not run as expected! Restart files are missing for some members!')
else:
    print('Nice. Spinup complete. Proceed to FREE!')

# if a postprocessing directory does not exist for this case, make one
output_path = '/glade/work/'+user+'/Projects/fractional-comp/data/processed_output/ensemble/'+case_name+'/'

if os.path.exists(output_path) == False:
    os.makedirs(output_path)

files = sorted(glob.glob('/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case_name+'/mem*/history/*.nc'))
DS = []
for file in files:
    ds = xr.open_dataset(file).isel({'ni':2}).drop(['ntrcr','ni','trcr','trcrn'])
    DS.append(ds)

ens_ds = xr.concat(DS, dim='member')
ens_ds = ens_ds.resample(time = '1D').mean()
ens_ds['hi'] = ens_ds.vice/ens_ds.aice

ens_ds.to_netcdf(output_path+'/postprocessed_ensemble.nc')
ens_ds.mean(dim='member').to_netcdf(output_path+'/postprocessed_ensemble_mean.nc')
ens_ds.std(dim='member', ddof=1).to_netcdf(output_path+'/postprocessed_ensemble_std.nc')

# check if files have been written
if len(glob.glob(output_path+'/postprocessed_ensemble*.nc')) == 3:
    print('Postprocessing complete! PROCESS FINISHED.')
else:
    print('Postprocessing failed!')