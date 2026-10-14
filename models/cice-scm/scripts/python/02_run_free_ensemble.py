#################################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import sys
import xarray as xr
import numpy as np

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the free case name and the spinup it is based on
case_name = sys.argv[1]
user = sys.argv[2]
spinup_case = sys.argv[3]

location = spinup_case.split('_')[1]
year_init = 2011
simulation_length = 1
ensemble_size = 30

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
    print('Model did not build correctly! Please rebuild model.')
else:
    print('Model has been built!')

#-----------------------------------------#
# 4. Run the case
#-----------------------------------------#
mem = 1
while mem <= ensemble_size:
    inst_string ='{0}'.format('%04d' % mem) 
    print('Running member '+inst_string+'...')

    #create history and restart directories for the run
    os.chdir(storage_dir)
    os.makedirs('mem' + inst_string + '/history/')
    os.makedirs('mem' + inst_string + '/restart/')

    # link the model executable for each member to the main one built for the case
    os.symlink(storage_dir+'/icepack','mem' + inst_string+'/icepack')
    
    # set the restart file used from spinup for this member
    spinup_dir = scratch_dir + '/ICEPACK_RUNS/' + spinup_case 
    restart_file = spinup_dir + '/mem'+inst_string+'/restart/iced.'+str(year_init)+'-01-01-00000.nc'
    
    # go the member directory and begin working with the namelist settings
    os.chdir(storage_dir+'/mem' + inst_string)
    namelist = f90nml.read(namelist_template)
   
    # set case settings 
    namelist['setup_nml']['year_init'] = year_init
    namelist['setup_nml']['npt'] = simulation_length * 8760
    namelist['setup_nml']['restart'] = True
    namelist['setup_nml']['runtype_startup'] = True
    namelist['setup_nml']['ice_ic'] = restart_file
    namelist['setup_nml']['dumpfreq'] = 'd'
    namelist['setup_nml']['restart_format'] = 'nc'
    namelist['setup_nml']['history_format'] = 'nc'
    if location != 'default':
        namelist['setup_nml']['input_lat'] = locs[location][0]
        namelist['setup_nml']['input_lon'] = locs[location][1]

    # change namelist parameters of note
    namelist['shortwave_nml']['r_snw'] = R_snw
    namelist['shortwave_nml']['dt_mlt'] = dt_mlt

    # set namelist forcing options
    namelist['forcing_nml']['data_dir'] = forcing_dir +location + '/free/'
    namelist['forcing_nml']['atm_data_type'] = 'JRA55'
    namelist['forcing_nml']['ocn_data_type'] = 'ISPOL'
    namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_' + member_string + '.txt'
    namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_' + member_string + '.txt'
    namelist['forcing_nml']['default_season'] = 'winter'
    namelist['forcing_nml']['lateral_flux_type'] = 'uniform_ice'
    namelist['forcing_nml']['update_ocn_f'] = True
    namelist['forcing_nml']['fyear_init'] = 2011
    # write namelist to needed file
    namelist.write('icepack_in',force=True)

    # run the member instance and dump output 
    comd = './icepack > icepack.out'
    os.system(comd)
    
    # advance to the next ensemble member 
    mem += 1

final_year = '{0}'.format('%04d' % (year_init + simulation_length))
check_restarts = glob.glob(storage_dir + '/mem*/restart/iced.'+final_year+'-01-01-00000.nc')

if len(check_restarts) != ensemble_size:
    print('HEY!! FREE run did not run as expected! Restart files are missing for some members!')
else:
    print('Nice. FREE complete. Proceed to postprocessing!')

# if a postprocessing directory does not exist for this case, make one
output_path = '/glade/work/'+user+'/Projects/fractional-comp/data/processed_output/free/'+case_name+'/'

if os.path.exists(output_path) == False:
    os.makedirs(output_path)

files = sorted(glob.glob('/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case_name+'/mem*/history/*.nc'))
DS = []
for file in files:
    ds = xr.open_dataset(file).isel({'ni':2}).drop_vars(['ntrcr','ni','trcr','trcrn'])
    DS.append(ds)

ens_ds = xr.concat(DS, dim='member')

ens_ds = ens_ds.resample(time = '1D').mean()
ens_ds['hi'] = ens_ds.vice/ens_ds.aice

ens_ds.to_netcdf(output_path+'/postprocessed_ensemble.nc')
ens_ds.mean(dim='member').to_netcdf(output_path+'/postprocessed_ensemble_mean.nc')
ens_ds.std(dim='member', ddof=1).to_netcdf(output_path+'/postprocessed_ensemble_std.nc')

# check if files have been written
if len(glob.glob(output_path+'/postprocessed_ensemble*.nc')) == 3:
    print('Postprocessing complete! PROCESSED FINISHED.')
else:
    print('Postprocessing failed!')