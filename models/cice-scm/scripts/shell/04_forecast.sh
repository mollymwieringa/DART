#!/bin/bash -l
### Job Name
#PBS -N forecast_Ba_DEC
### Charging account
#PBS -A P93300065
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=4GB
### Allow job to run up to 30 minutes
#PBS -l walltime=01:00:00 
#PBS -l job_priority=economy
### Route the job to the casper queue
#PBS -q main 
### Join output and error streams into single file
#PBS -j oe
#PBS -m ae
#PBS -M mmw906@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
conda activate cice-scm-da

python ../python/04b_forecast.py SIT_Barents_DEC_BNRH_init mollyw SIT_Barents_BNRH 2011 12 15 4
python ../python/04b_forecast.py SIC_Barents_DEC_BNRH_init mollyw SIC_Barents_BNRH 2011 12 15 4
python ../python/04b_forecast.py FB_Barents_DEC_BNRH_init mollyw FB_Barents_BNRH 2011 12 15 4