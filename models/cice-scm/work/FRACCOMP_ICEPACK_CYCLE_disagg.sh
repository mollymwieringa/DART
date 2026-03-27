#!/bin/bash -l
### Job Name
#PBS -N SibChuk_disaggregation_cycling
### Charging account
#PBS -A P54048000
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
#PBS -M mollyw@ucar.edu

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load Python module and activate NPL environment
conda activate icepack

### Run analysis script
python FRACCOMP_ICEPACK_CYCLE.py SibChuk SIC disaggregation 3 2011 1 2 2011 12 31