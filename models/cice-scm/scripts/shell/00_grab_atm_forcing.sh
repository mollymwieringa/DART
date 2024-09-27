#!/bin/bash -l
### Job Name
#PBS -N grab_atm_forcings
### Charging account
#PBS -A UWAS0083
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=100GB
### Allow job to run up to 30 minutes
#PBS -l walltime=1:00:00
#PBS -l job_priority=economy
### Route the job to the economy queue
#PBS -q main
### Join output and error streams into single file
#PBS -j oe
### send emails on abort and exit
#PBS -m ae
#PBS -M mmw906@uw.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
# module load ncarenv python
conda activate cice-scm-da

### Run analysis script
python ../python/00_grab_atm_forcing.py 2000 2000 spinup 2