
import sys 
import numpy as np
import os
import xarray as xr

# Assign command line arguments
loc = str(sys.argv[1]) # SibChuk
obs_type = str(sys.argv[2]) # SIC
state_type = str(sys.argv[3])
g = float(sys.argv[4])# 0.1


# Handle file system reading / creating 
base_path = '/glade/u/home/mollyw/work/Projects/'
file = base_path + 'cice-scm-da/data/experiment_output/regime_assimilation/free/free_'+loc+'/postprocessed_ensemble.nc'
ds = xr.open_dataset(file).isel(time=slice(-366, -1))

project_path = base_path + 'fractional-comp/experiment_output/offline/' #fixed_error_01/fix_CA_SIT_aicen/'
case_path = project_path + loc + '/' + obs_type + '_' + state_type
os.makedirs(case_path, exist_ok=True)

# Gather relevant loop / sampling variables
day_count = 0
rand_int = np.random.randint(0,29)
# rand_int = 20
for t in ds.time:
    day_count += 1
    day_path = case_path + '/day_{:03d}'.format(day_count)
    os.makedirs(day_path, exist_ok=True)

    # Compute the aggregate values based on the observations type. 
    # Assume that aggregates exist only where there is at least 15% sea ice concentration 
    sic = np.sum(ds.sel(time=t).aicen.values, axis=1)
    sic = np.where(sic > 1e-6, sic, 0)

    if obs_type == 'SIC':
        aggregates = np.sum(ds.sel(time=t).aicen.values, axis=1)
    elif obs_type == 'SIT':
        div = np.where(sic > 0, sic, -9999.0)
        aggregates = np.sum(ds.sel(time=t).vicen.values, axis=1) / div
        aggregates = np.where(aggregates >= 0.0, aggregates, 0.0)
    elif obs_type == 'SND':
        div = np.where(sic > 0, sic, -9999.0)
        aggregates = np.sum(ds.sel(time=t).vsnon.values, axis=1) / div
        aggregates = np.where(aggregates >= 0.0, aggregates, 0.0)
    else:
        print('obs type not recognized, exiting')
        sys.exit()
    
    # Isolate the state variable that is being updated by the relevant observation type
    ds_day = ds.sel(time=t)[state_type].values

    # Extract observation value at a random ensemble member and add noise
    y = aggregates[rand_int]
    noise = np.random.normal(0, g, 1)
    obs = y + noise[0]

    # check that the observation is within physical bounds
    if (obs_type == 'SIC') & (obs > 1.0):
        obs = 1.0
    elif obs < 0.0:
        obs = 0.0
    else:
        pass

    # Save the prior and observation info to text files for assimilation 
    np.savetxt(f'obs_info.txt', [obs, g, rand_int, y])
    np.savetxt(f'true_cats_info.txt', ds_day[rand_int])

    ds_day = np.delete(ds_day, rand_int, axis=0)
    aggregates = np.delete(aggregates, rand_int, axis=0)
    np.savetxt(f'prior_ensemble.txt', ds_day)
    np.savetxt(f'obs_prior_ensemble.txt', aggregates)

    # if all ensemble members in all categories are zero, skip this day
    if np.sum(ds_day) == 0.0:
        comd = f'echo "All members in all categories are zero, skipping day " {day_count}'
        os.system(comd)
    else:
        comd = f'echo "Assimilating day " {day_count}'
        os.system(comd)
        
        # Calculate the observation increments 
        comd = './obs_increments'
        os.system(comd)

        # if the obs_posterior file exists, run the state_regression
        if os.path.isfile('obs_post_ensemble_bnrhf.txt'):
            # Do the BNRHF regression
            comd = f'cp input_bnrhf.nml input.nml'
            os.system(comd)
            comd = './state_regression'
            os.system(comd)

            # move the BNRHF output files to the storage directory
            comd = f'mv obs_post_ensemble_bnrhf.txt '+ day_path +'/obs_post_ensemble_bnrhf.txt'
            os.system(comd)
            comd = f'mv ens_post_disaggregation_bnrh.txt '+ day_path +'/ens_post_disaggregation_bnrh.txt'
            os.system(comd)
            comd = f'mv ens_post_relativefrac_bnrh.txt '+ day_path +'/ens_post_relativefrac_bnrh.txt'
            os.system(comd)
            if (os.path.isfile('alpha_prior_bnrh.txt')):
                comd = f'mv alpha_prior_bnrh.txt ' + day_path +'/alpha_prior_bnrh.txt'
                os.system(comd)
            comd = f'mv alpha_posterior_bnrh.txt ' + day_path +'/alpha_posterior_bnrh.txt'
            os.system(comd)
            comd = f'mv ens_post_probit_raw_bnrh.txt '+ day_path +'/ens_post_probit_raw_bnrh.txt'
            os.system(comd)
            comd = f'mv ens_post_probit_postprocess_bnrh.txt '+ day_path +'/ens_post_probit_postprocess_bnrh.txt'
            os.system(comd)
            comd = f'mv ens_post_linear_raw_bnrh.txt '+ day_path +'/ens_post_linear_raw_bnrh.txt'
            os.system(comd)
            comd = f'mv ens_post_linear_postprocess_bnrh.txt '+ day_path +'/ens_post_linear_postprocess_bnrh.txt'
            os.system(comd)
            comd = f'mv obs_increments_bnrhf.txt '+ day_path +'/obs_increments_bnrhf.txt'
            os.system(comd)
        else:
            comd = f'echo Observation incrementing for BNRHF did not succeed at time '+str(t.values)+', therefore I cannot update the state!'
            os.system(comd)
            # print('Observation incrementing for BNRHF did not succeed at time '+str(t.values)+', therefore I cannot update the state!')

        if os.path.isfile('obs_post_ensemble_kde.txt'):
            # Do the KDE regression
            comd = f'cp input_kqcef.nml input.nml'
            os.system(comd)
            comd = './state_regression'
            os.system(comd)

            # move the KDE output files to the storage directory 
            comd = f'mv obs_post_ensemble_kde.txt '+ day_path +'/obs_post_ensemble_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_disaggregation_kde.txt '+ day_path +'/ens_post_disaggregation_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_relativefrac_kde.txt '+ day_path +'/ens_post_relativefrac_kde.txt'
            os.system(comd)
            if (os.path.isfile('alpha_prior_kde.txt')):
                comd = f'mv alpha_prior_kde.txt ' + day_path +'/alpha_prior_kde.txt'
                os.system(comd)
            comd = f'mv alpha_posterior_kde.txt ' + day_path +'/alpha_posterior_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_probit_raw_kde.txt '+ day_path +'/ens_post_probit_raw_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_probit_postprocess_kde.txt '+ day_path +'/ens_post_probit_postprocess_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_linear_raw_kde.txt '+ day_path +'/ens_post_linear_raw_kde.txt'
            os.system(comd)
            comd = f'mv ens_post_linear_postprocess_kde.txt '+ day_path +'/ens_post_linear_postprocess_kde.txt'
            os.system(comd)
            comd = f'mv obs_increments_kde.txt '+ day_path +'/obs_increments_kde.txt'
            os.system(comd)
        else:
            comd = f'echo Observation incrementing for KQCEF did not succeed at time '+str(t.values)+', therefore I cannot update the state!'
            os.system(comd)
            # print('Observation increments for KQCEF did not succeed at time '+str(t.values)+', therefore I cannot update the state!')

    # move the observation / verification info to the storage directory
    comd = f'mv obs_info.txt '+ day_path +'/obs_info.txt'
    os.system(comd)
    comd = f'mv true_cats_info.txt '+ day_path +'/true_cats_info.txt'
    os.system(comd)

    # Once assimilation is complete, move the prior files to the storage directory 
    comd = f'mv prior_ensemble.txt '+ day_path +'/prior_ensemble.txt'
    os.system(comd)
    comd = f'mv obs_prior_ensemble.txt '+ day_path +'/obs_prior_ensemble.txt'
    os.system(comd)

print('done!')