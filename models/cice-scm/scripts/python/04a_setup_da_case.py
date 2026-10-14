###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import xarray as xr
import numpy as np
import sys

from datetime import datetime

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the free cases name and the spinup it is based on
case_name = sys.argv[1] # SIT_f1_NORM_test
user = sys.argv[2]
spinup_case = sys.argv[3] # spinup_test

ensemble_size = 30
year_init = 2011
first_assim_time = datetime(2011, 1, 2)
model_init_time = datetime(2011, 1, 1)
location = spinup_case.split('_')[1]

# set the relevant directory names
# where is this project repo located on your system?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# where is the Icepack installation located on your system?
icepack_dir = '/glade/work/'+user+'/Icepack/'
# where is the DART installation located on your system?
dart_dir = '/glade/work/'+user+'/dart_manhattan/'
# where is the scratch directory Icepack will write output to?
scratch_dir = '/glade/derecho/scratch/'+user+'/'
# where do your forcing files live? should have spinup, free subdirectories
forcing_dir = project_dir + '/data/forcings/'
# where does your Icepack output get written to? 
storage_dir = scratch_dir + '/ICEPACK_RUNS/' + case_name
# where is your icepack_in template file?
namelist_template = icepack_dir + case_name + '/icepack_in'

# set the machine name and compiler
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

days_to_assim = first_assim_time - model_init_time
assim_date_str = '{:04d}-{:02d}-{:02d}'.format(first_assim_time.year, 
                                               first_assim_time.month, 
                                               first_assim_time.day)
simulation_length = 24 * days_to_assim.days
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
    # AssertionError('Model did not build correctly! Please rebuild model.')
else:
    print('Model with following perturbations has been built:', perturb)

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
    restart_file = spinup_dir + '/mem' + member_string + '/restart/iced.'+str(year_init)+'-01-01-00000.nc'

    # begin working on an individual ensemble member
    os.chdir(storage_dir+'/mem' + inst_string)
    namelist = f90nml.read(namelist_template)
    

    # set case settings 
    namelist['setup_nml']['year_init'] = year_init
    namelist['setup_nml']['npt'] = simulation_length
    namelist['setup_nml']['restart'] = True
    namelist['setup_nml']['runtype_startup'] = True
    namelist['setup_nml']['ice_ic'] = restart_file    
    namelist['setup_nml']['dumpfreq'] = 'd'
    namelist['setup_nml']['restart_format'] = 'nc'
    namelist['setup_nml']['history_format'] = 'nc'
    namelist['setup_nml']['conserv_check'] = True
    if location != 'default':
        namelist['setup_nml']['input_lat'] = locs[location][0]
        namelist['setup_nml']['input_lon'] = locs[location][1]

    # set parameter values for this member 
    namelist['shortwave_nml']['r_snw'] = R_snw
    namelist['shortwave_nml']['dt_mlt'] = dt_mlt

    # set forcing directory for this member
    namelist['forcing_nml']['data_dir'] = forcing_dir + location + '/free/'     
    namelist['forcing_nml']['atm_data_type'] = 'JRA55'
    namelist['forcing_nml']['ocn_data_type'] = 'ISPOL'
    namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_' + member_string + '.txt'
    namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_' + member_string + '.txt'
    namelist['forcing_nml']['default_season'] = 'winter'
    namelist['forcing_nml']['lateral_flux_type'] = 'uniform_ice'
    namelist['forcing_nml']['update_ocn_f'] = True
    namelist['forcing_nml']['fyear_init'] = year_init
    # write namelist to needed file
    namelist.write('icepack_in',force=True)

    # run the member instance and dump output 
    comd = './icepack > icepack.out'
    os.system(comd)
    
    # advance to the next ensemble member 
    mem += 1

# collect a list of all the final (year 2011) restart files
check_restarts = glob.glob(storage_dir + '/mem*/restart/iced.'+str(year_init)+'-01-02-00000.nc')

if len(check_restarts) != ensemble_size:
    print('HEY!! Assimilation run did not run as expected! Restart files are missing for some members!')
else:
    print('Nice. Assimilation start is complete. Proceed to cycling!')
