import glob
import os
import xarray as xr
import sys

case_name = sys.argv[1]
user = sys.argv[2]

output_path = '/glade/work/'+user+'/Projects/cice-scm-da/data/processed/ensemble/'+case_name+'/'
if os.path.exists(output_path) == False:
    os.mkdir(output_path)

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
if len(glob.glob(output_path+'/postprocessed_ens*.nc')) == 3:
    print('Postprocessing complete! PROCESSED FINISHED.')
else:
    print('Postprocessing failed!')