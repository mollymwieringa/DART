###############################################################
# IMPORT NECESSARY LIBRARIES                                  #
###############################################################
import sys
import os
import shutil
import glob
import f90nml
import numpy as np
import xarray as xr
import pandas as pd
from time import sleep, perf_counter
from datetime import datetime, timedelta
###############################################################
# SET MANUAL INPUTS HERE                                      #
###############################################################
loc = sys.argv[1] #must be CentArc Barents SibChuk
obs_type=sys.argv[2]
regression_type = sys.argv[3]
truth_member = sys.argv[4]

# set the assimilation dates
first_assim_date = datetime(int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]))
end_assim_date = datetime(int(sys.argv[8]),int(sys.argv[9]),int(sys.argv[10]))

ens_size = 30

case = loc+'_'+obs_type+'_'+regression_type
spinup_case = 'SPINUP_SibChuk'
###############################################################
# DEFINE HELPER FUNCTIONS                                     #
###############################################################

# output files are 
#   - all the txt files for each day (include all the ensemble members in them)
#   - restart files for the next timesteps (for each ensemble member and each day)
# Required directories for Icepack (OUTPUT ICEPACK FILES)
#   - Icepack case directory (scratch/ICEPACK_RUNS/case/mem_string)
#   - history and restart files for each member are in this case directory
# Specialized directories (OUTPUT TEXT FILES)
#   - work/Projects/fractional-comp/experiment_output/online/
#   - subdirectories for each regime >> each obs/state combo >> each date


def run_assimilation(icepack_path, storage_path, obs_type, regression_type, assim_date, ens_size, truth_member):

    # set up need info about this assimilation cycle
    date_str = '{0}'.format('%04d'%assim_date.year) + '-{0}'.format('%02d'%assim_date.month) + '-{0}'.format('%02d'%assim_date.day)

    # establish paths
    day_path = storage_path + '/'+case+'/'+date_str+'/'
    os.makedirs(day_path, exist_ok=True)

    # bring in DART namelist file
    comd = 'cp /glade/work/mollyw/dart_manhattan/models/cice-scm/scripts/templates/DART_input_cycle.nml input.nml'
    os.system(comd)
    dart_nml = f90nml.read('input.nml')
    dart_nml['obs_increments_nml']['bounded_above'] = True
    dart_nml['obs_increments_nml']['upper_bound'] = 1.0
    dart_nml['state_regression_nml']['obs_bounded_above'] = True
    dart_nml['state_regression_nml']['obs_upper_bound'] = 1.0
    dart_nml['state_regression_nml']['regression_type'] = regression_type
    dart_nml.write('input.nml', force=True)


    # bring in observation info
    comd = 'cp '+storage_path+'/observations/'+obs_type+'/obs_info_'+date_str+'.txt obs_info.txt'
    os.system(comd)
    comd = 'cp '+storage_path+'/observations/'+obs_type+'/true_cats_info_'+date_str+'.txt true_cats_info.txt'
    os.system(comd)
    comd = 'cp '+storage_path+'/'+case+'/ensemble/prior_ensemble_'+date_str+'.txt prior_ensemble.txt'
    os.system(comd)
    comd = 'cp '+storage_path+'/'+case+'/ensemble/obs_prior_ensemble_'+date_str+'.txt obs_prior_ensemble.txt'
    os.system(comd)

    # run assimilation steps
    comd = './obs_increments > obs_increment_output'
    os.system(comd)

    if (os.path.isfile('obs_post_ensemble_bnrhf.txt') is False):
        print('Obs. incrementing did not finish correctly! Cannot perform state regression. Exiting...')
        output = 0
        sys.exit()
    elif (abs(np.mean(np.loadtxt('obs_increments_bnrhf.txt')) - 0.0 ) < 1e-12):
        print('Obs. incrementing on '+date_str+' did not adjust the SIC ensemble at all! Skipping state regression and moving on to forecast step...')
        output = 1
    else:
        comd = './state_regression > state_regression_output'
        os.system(comd)
        print(regression_type)
        # update the restart file for icepack and tidy up
        if (os.path.isfile('ens_post_'+regression_type+'_bnrh.txt') is False):
            print('State regression did not finish correctly. Exiting...')
            output = 0
            sys.exit()
        else:
            # write all the ensemble members output to restart files 
            ens_post = np.loadtxt('ens_post_'+regression_type+'_bnrh.txt')
            assert np.shape(ens_post) == (ens_size-1,5), 'error loading ens. posterior data! expecting dimensions of '+str((ens_size-1,5))+'! got dimensions of '+str(np.shape(ens_post))+'.'
            mem_counter = 0
            for mem in range(ens_size):
                if mem+1 == int(truth_member):
                    print('Skipping truth member!')
                else:
                    # identify the restart file for this date and time
                    inst_string = '{0}'.format('%04d'%(mem+1))
                    restart_file = icepack_path+'/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc'

                    # make a copy of this restart file
                    comd = f'cp '+restart_file+' '+icepack_path+'/mem'+inst_string+'/restart/iced.'+date_str+'-00000_original.nc'
                    os.system(comd)

                    # replace the aicen categories in this restart file
                    restart_ds = xr.load_dataset(restart_file)
                    aicen_new = np.array(ens_post[mem_counter,:])
                    restart_ds['aicen'][:,2] = aicen_new
                    restart_ds.to_netcdf(restart_file)
                    mem_counter += 1

            # move state_regression output to project output directory 
            comd=f'mv obs_info.txt '+day_path+'/obs_info.txt'
            os.system(comd)
            comd=f'mv obs_prior_ensemble.txt '+day_path+'/obs_prior_ensemble.txt'
            os.system(comd)
            comd=f'mv obs_post_ensemble_bnrhf.txt '+day_path+'/obs_post_ensemble.txt'
            os.system(comd)
            comd = f'mv obs_increments_bnrhf.txt '+day_path+'/obs_increments.txt'
            os.system(comd)
            comd= f'mv true_cats_info.txt '+day_path+'/true_cats_info.txt'
            os.system(comd)
            comd=f'mv prior_ensemble.txt '+day_path+'/prior_ensemble.txt'
            os.system(comd)
            comd=f'mv ens_post_'+regression_type+'_bnrh.txt '+day_path+'/ens_post_'+regression_type+'.txt'
            os.system(comd)
            comd=f'mv state_regression_output '+day_path+'/state_regression_output'
            os.system(comd)
            comd=f'mv obs_increment_output '+day_path+'/obs_increment_output'
            os.system(comd)

            if regression_type in ['probit_postprocess', 'linear_postprocess']:
                comd=f'mv ens_post_'+regression_type[:-12]+'_raw_bnrh.txt '+day_path+'/ens_post_'+regression_type[:-12]+'_raw.txt'
                os.system(comd)
            if os.path.isfile('alpha_prior_bnrh.txt'):
                comd=f'mv alpha_prior_bnrh.txt '+day_path+'/alpha_prior.txt'
                os.system(comd)
            if os.path.isfile('alpha_posterior_bnrh.txt'):
                comd=f'mv alpha_posterior_bnrh.txt '+day_path+'/alpha_posterior.txt'
                os.system(comd)
            output = 1

    return output

def run_icepack(icepack_path, storage_path, obs_type, assim_date, truth_member):
    
    # parse time information
    date_str = '{0}'.format('%04d'%assim_date.year) + '-{0}'.format('%02d'%assim_date.month) + '-{0}'.format('%02d'%assim_date.day)
    next_day = assim_date + timedelta(days=1)
    next_day_str = '{0}'.format('%04d'%next_day.year) + '-{0}'.format('%02d'%next_day.month) + '-{0}'.format('%02d'%next_day.day) 

    mem = 1
    while mem <= ens_size:
        inst_string = '{0}'.format('%04d'%mem) 
        os.chdir(icepack_path + '/mem'+inst_string+'/')

        # update namelist with date info
        icepack_nml = f90nml.read('icepack_in')
        icepack_nml['setup_nml']['ice_ic'] = icepack_path + '/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc'
        icepack_nml['setup_nml']['runtype_startup'] = False
        icepack_nml.write('icepack_in', force=True)

        # run the model
        comd = './icepack > icepack.out'
        os.system(comd)
        sleep(2)

        # check if the simulation finished correctly
        restart_file = 'restart/iced.'+next_day_str+'-00000.nc'
        if os.path.exists(restart_file) is False:
            print('Icepack did not create necessary restart file! Process stopping.')
            output = 0
            sys.exit()
        else:
            output = 1
        
        mem += 1

    # write model forecast to next obs_prior_ensemble_DATE.txt and prior_ensemble_DATE.txt for this member to the storage directory
    with open('new_prior_ensemble.txt', 'a') as f_state:
        mem = 1
        while mem <= ens_size:
            if mem == int(truth_member):
                print('Skipping truth member!')
            else:
                inst_string = '{0}'.format('%04d'%mem) 
                restart_file = icepack_path + '/mem'+inst_string+'/restart/iced.'+next_day_str+'-00000.nc'
                restart_ds = xr.load_dataset(restart_file).isel({'ni':2})
                prior_values = restart_ds.aicen.values
                for i, val in enumerate(prior_values):
                    if i < 4:
                        f_state.write(str(val)+'   ')
                    else:
                        f_state.write(str(val)+'\n')
            mem += 1

    with open('new_obs_prior_ensemble.txt', 'a') as f_obs:
        mem = 1
        while mem <= ens_size:
            if mem == int(truth_member):
                print('Skipping truth member!')
            else:
                inst_string = '{0}'.format('%04d'%mem) 
                restart_file = icepack_path + '/mem'+inst_string+'/restart/iced.'+next_day_str+'-00000.nc'
                restart_ds = xr.load_dataset(restart_file).isel({'ni':2})
                if obs_type == 'SIC':
                    obs_prior_value = restart_ds.aicen.sum(dim='ncat').values
                elif obs_type == 'SIT':
                    obs_prior_value = restart_ds.vicen.sum(dim='ncat').values/restart_ds.aicen.sum(dim='ncat').values
                else:
                    print('obs type not recognized!')
                f_obs.write(str(obs_prior_value)+'\n')
            mem += 1

        # second, clean up 
        # shutil.copy('history/icepack.h.{0}{1}{2}.nc'.format('%04d'%assim_date.year, '%02d'%assim_date.month, '%02d'%assim_date.day), icepack_path + '/mem'+inst_string+'/forecasts/'+date_str+'/icepack.h.'+date_str+'_'+inst_string+'.nc')
        # shutil.move('ice_diag.full_ITD', icepack_path + '/mem'+inst_string+'/forecasts/'+date_str+'/ice_diag.full_ITD_'+inst_string)
        # os.remove('icepack.out')
        # comd = 'rm ice_diag.*'
        # os.system(comd)

    # move the new ensemble text files
    comd = 'mv new_prior_ensemble.txt '+storage_path+'/'+case+'/ensemble/prior_ensemble_'+next_day_str+'.txt'
    os.system(comd)
    comd = 'mv new_obs_prior_ensemble.txt '+storage_path+'/'+case+'/ensemble/obs_prior_ensemble_'+next_day_str+'.txt'
    os.system(comd)

    return output

def cycle(case, assim_date, ens_size, obs_type, truth_member):

    icepack_path =  '/glade/derecho/scratch/mollyw/ICEPACK_RUNS/'+case+'/'
    assim_path = icepack_path + '/assim_dir/'
    dart_path = '/glade/work/mollyw/dart_manhattan//models/cice-scm/work/'
    storage_path = '/glade/work/mollyw/Projects/fractional-comp/data/experiment_output/online/'
    
    regression_type = case[12:]

    # go to assim directory and run assimilation
    #if assim_path does not exist, make it
    if not os.path.exists(assim_path):
        os.makedirs(assim_path)
    os.chdir(assim_path)

    # copy over executables from the dart path
    comd = 'cp '+dart_path+'obs_increments .'
    os.system(comd)
    comd = 'cp '+dart_path+'state_regression .'
    os.system(comd)
    run_assimilation(icepack_path, storage_path, obs_type, regression_type, assim_date, ens_size, truth_member)

    # go to Icepack and run forecast for this timestep
    os.chdir(icepack_path)
    run_icepack(icepack_path, storage_path, obs_type, assim_date, truth_member)

    return

def setup(case, spinup_case, obs_type, ens_size):
    """Goal with this is to set up the bases of the experimental case. Assumes that a free run 
    and associated observations have already been created"""

    os.chdir('/glade/work/mollyw/Icepack/')
    comd = './icepack.setup -c '+case+' -m derecho -e inteloneapi'
    os.system(comd)

    # go to the case directory
    os.chdir(case)

    # build the case
    comd = './icepack.build'
    os.system(comd)

    # check that the case was built correctly
    storage_dir = '/glade/derecho/scratch/mollyw/ICEPACK_RUNS/' + case
    if os.path.exists(storage_dir) is False:
        print('Model did not build correctly! Please rebuild model.')

    if os.path.exists('/glade/work/mollyw/Projects/fractional-comp/data/experiment_output/online/'+case+'/ensemble/') is False:
        os.makedirs('/glade/work/mollyw/Projects/fractional-comp/data/experiment_output/online/'+case+'/ensemble/')

    mem = 1
    while mem <= ens_size:
        inst_string ='{0}'.format('%04d' % mem) 
        print('Running member '+inst_string+'...')
        restart_file = '/glade/derecho/scratch/mollyw/ICEPACK_RUNS/'+spinup_case+'/mem'+inst_string+'/restart/iced.2011-01-01-00000.nc'
        ds = xr.load_dataset(restart_file)
        ds.attrs['istep1'] = 0.0
        ds.attrs['time'] = 0.0
        ds.to_netcdf(restart_file)

        # create history and restart directories for the run
        os.chdir(storage_dir)
        os.makedirs('mem' + inst_string + '/history/')
        os.makedirs('mem' + inst_string + '/restart/')

        # link the model executable for each member to the main one built for the case
        os.symlink(storage_dir+'/icepack','mem' + inst_string+'/icepack')

        # begin working on an individual ensemble member
        os.chdir(storage_dir+'/mem' + inst_string)

        # copy in the relevant namelist
        # read namelist template
        namelist = f90nml.read('/glade/work/mollyw/Projects/fractional-comp/data/templates/icepack_in.setup')
        namelist['setup_nml']['ice_ic'] = restart_file
        namelist['setup_nml']['input_lat'] = 1.3183906501
        namelist['setup_nml']['input_lon'] = 3.0446500856
        namelist['setup_nml']['runtype_startup'] = True
        namelist['forcing_nml']['data_dir'] = '/glade/work/mollyw/Projects/cice-scm-da/data/forcings/SibChuk/free/'
        namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_'+inst_string+'.txt'
        namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_'+inst_string+'.txt'
        namelist.write('icepack_in',force=True)
        
        comd = './icepack > icepack.out'
        os.system(comd)

        mem += 1
    
    # write model forecast to next obs_prior_ensemble_DATE.txt and prior_ensemble_DATE.txt for this member to the storage directory
    os.chdir(storage_dir)
    with open('first_prior_ensemble.txt', 'a') as f_state:
        mem = 1
        while mem <= ens_size:
            if mem == int(truth_member):
                print('Skipping truth member!')
            else:
                inst_string ='{0}'.format('%04d' % mem) 
                new_restart_file = '/glade/derecho/scratch/mollyw/ICEPACK_RUNS/'+case+'/mem'+inst_string+'/restart/iced.2011-01-02-00000.nc'
                restart_ds = xr.load_dataset(new_restart_file).isel({'ni':2})
                prior_values = restart_ds.aicen.values
                for i, val in enumerate(prior_values):
                    if i < 4:
                        f_state.write(str(val)+'   ')
                    else:
                        f_state.write(str(val)+'\n')
                restart_ds.close()
            mem += 1
    
    with open('first_obs_prior_ensemble.txt', 'a') as f_obs:
        mem = 1
        while mem <= ens_size:
            if mem == int(truth_member):
                print('Skipping truth member!')
            else:
                inst_string ='{0}'.format('%04d' % mem) 
                new_restart_file = '/glade/derecho/scratch/mollyw/ICEPACK_RUNS/'+case+'/mem'+inst_string+'/restart/iced.2011-01-02-00000.nc'
                restart_ds = xr.load_dataset(new_restart_file).isel({'ni':2})
                if obs_type == 'SIC':
                    obs_prior_value = restart_ds.aicen.sum(dim='ncat').values
                elif obs_type =='SIT':
                    obs_prior_value = restart_ds.vicen.sum(dim='ncat').values/restart_ds.aicen.sum(dim='ncat').values
                else:
                    print('obs_type not recognized!')
                f_obs.write(str(obs_prior_value)+'\n')
                restart_ds.close()
            mem += 1

    # move the new ensemble text files
    comd = 'mv '+storage_dir+'/first_prior_ensemble.txt /glade/work/mollyw/Projects/fractional-comp/data/experiment_output/online/'+case+'/ensemble/prior_ensemble_2011-01-02.txt'
    os.system(comd)
    comd = 'mv '+storage_dir+'/first_obs_prior_ensemble.txt /glade/work/mollyw/Projects/fractional-comp/data/experiment_output/online/'+case+'/ensemble/obs_prior_ensemble_2011-01-02.txt'
    os.system(comd)

    # write the first round of prior ensemble files to /glade/work/mollyw/Projects/fractional-comp/experiment_output/online/'+case+'/ensemble/'

    return

###############################################################
# PERFORM CYCLING                                             #
###############################################################

start = perf_counter()
setup(case, spinup_case, obs_type, ens_size)
end = perf_counter()
print(f'Setup for the case took {(end-start)/60:0.4f} minutes.')

start = perf_counter()
datelist = pd.date_range(start=first_assim_date, end=end_assim_date).to_pydatetime()

for assim_date in datelist:
    cycle(case, assim_date, ens_size, obs_type, truth_member)
end = perf_counter()
print(f'Execution for cycling took {(end - start)/60:0.4f} minutes.')


# 1. Have truth data for each cycling day from the ensemble
# 2. Generate observations from the truth member 
# 3. For each day of the year:
#   - get the prior from the restart file for that day
#   - write the prior to text files (to be read by obs_increment)
#   - use obs_increment to calculate the obs_space adjustment
#   - use the obs_increment text output to run state_regression
#   - write the posterior in text files to icepack restart file
#   - move all text file output to daily folder for this case
#   - run a forecast with Icepack to the next timestep
