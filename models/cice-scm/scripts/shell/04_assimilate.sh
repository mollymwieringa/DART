#!/bin/bash -l
### Job Name
#PBS -N run_assimilation
### Charging account
#PBS -A P93300065
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=4GB
### Allow job to run up to 30 minutes
#PBS -l walltime=11:00:00 
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

### Run analysis script
python ../python/04a_setup_da_case.py SIC_SibChuk_relativefrac_3 mollyw SPINUP_SibChuk
python ../python/04b_cycle.py SIC_SibChuk_relativefrac_3 mollyw FREE_SibChuk 3 cice 2011 1 2 2011 12 31 SAT_SEAICE_AGREG_CONCENTR
python ../python/05_postprocess.py SIC_SibChuk_relativefrac_3 mollyw 13 all
