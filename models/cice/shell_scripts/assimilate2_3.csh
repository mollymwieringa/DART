#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN ICE_ASSIMILATE"

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

set PARAMETER_ESTIMATION = "FALSE"
set use_mean = "FALSE"
set cice_parameters = ( r_snw Cf )
# these aliases look like a module bug workaround.  if you get
# errors running the nco commands you can try to uncomment these, 
# but they should NOT be necessary.  feel free to remove these 
# lines if you test this script and it just works.
#set ncodir = "/glade/u/apps/ch/opt/nco/4.6.2/gnu/6.3.0/bin/"
#alias ncks $ncodir/ncks
#alias ncwa $ncodir/ncwa
#alias ncrename $ncodir/ncrename
#alias ncap2 $ncodir/ncap2
module load nco

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case ys*:
         # NCAR "yellowstone"
         set      MOVE = '/bin/mv -v'
         set      COPY = '/bin/cp -v --preserve=timestamps'
         set      LINK = '/bin/ln -vs'
         set    REMOVE = '/bin/rm -rf'
         set LAUNCHCMD = mpirun.lsf
         set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
         setenv MP_DEBUG_NOTIMEOUT yes
      breaksw
   case r*:
     # cheyenne has nodes like r1i0n14
         set      MOVE = '/usr/bin/mv -v'
         set      COPY = '/usr/bin/cp -v --preserve=timestamps'
         set      LINK = '/usr/bin/ln -vs'
         set    REMOVE = '/usr/bin/rm -rf'
         set LAUNCHCMD = 'mpiexec_mpt omplace'
   breaksw

   case linux_system_with_utils_in_other_dirs*:
         # example of pointing this script at a different set of basic commands
         set      MOVE = '/usr/local/bin/mv -v'
         set      COPY = '/usr/local/bin/cp -v --preserve=timestamps'
         set      LINK = '/usr/local/bin/ln -vs'
         set    REMOVE = '/usr/local/bin/rm -fr'
         set LAUNCHCMD = mpirun.lsf
      breaksw
   default:
         # NERSC "hopper"
         set      MOVE = 'mv -v'
         set      COPY = 'cp -v --preserve=timestamps'
         set      LINK = 'ln -vs'
         set    REMOVE = 'rm -fr'
         set LAUNCHCMD = "aprun -n 1" #$NTASKS"
      breaksw
endsw

# The bogus strings get replaced when CESM_DART_config is run
setenv    CASEROOT BOGUSCASEROOT
setenv BASEOBSROOT BOGUSBASEOBSDIR
#-------------------------------------------------------------------------
# Get the case-specific variables
#-------------------------------------------------------------------------

cd ${CASEROOT} || exit 1
setenv CASE          `./xmlquery --value CASE     `
setenv EXEROOT       `./xmlquery --value EXEROOT  `
setenv RUNDIR        `./xmlquery --value RUNDIR   `
setenv ensemble_size `./xmlquery --value NINST_ICE`
setenv ICE_COMPONENT `./xmlquery --value COMP_ICE`
setenv archive       `./xmlquery --value DOUT_S_ROOT`
setenv CONTINUE_RUN  `./xmlquery --value CONTINUE_RUN`
# Check to make sure we are running what we are supporting

if ( $ICE_COMPONENT != 'cice' ) then
   echo 'ERROR: This assimilate.csh file is specifically for "cice".'
   echo 'ERROR: the ice component for this case is "'${ICE_COMPONENT}'".'
   exit 1
endif

cd ${RUNDIR}

#-------------------------------------------------------------------------
# Determine time of model state ... from the last coupler restart file name
# of the form "./${CASE}.cpl_*.r.YYYY-MM-DD-SSSSS.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set REST_FILE = `ls -1 $CASE.cpl_*.r.*.nc | tail -n 1`
set HIST_FILE = `ls -1 $CASE.cice_*.h.*.nc | tail -n 1`
set modelFILE = $REST_FILE:r
set obsFILE = $HIST_FILE:r

set CPL_DATE_EXT = `echo $modelFILE:e`
set CPL_DATE     = `echo $modelFILE:e | sed -e "s#-# #g"`
set OBS_DATE_EXT = `echo $obsFILE:e`
set OBS_DATE     = `echo $obsFILE:e | sed -e "s#-# #g"`

set CPL_YEAR     = `echo $CPL_DATE[1] | bc`
set OBS_YEAR     = `echo $OBS_DATE[1] | bc`
set CPL_MONTH    = `echo $CPL_DATE[2] | bc`
set OBS_MONTH    = `echo $OBS_DATE[2] | bc`
set CPL_DAY      = `echo $CPL_DATE[3] | bc`
set OBS_DAY      = `echo $OBS_DATE[3] | bc`
set CPL_SECONDS  = `echo $CPL_DATE[4] | bc`
set CPL_HOUR     = `echo $CPL_DATE[4] / 3600 | bc`

echo "valid time of model is $CPL_YEAR $CPL_MONTH $CPL_DAY $CPL_SECONDS (seconds)"
echo "valid time of model is $CPL_YEAR $CPL_MONTH $CPL_DAY $CPL_HOUR (hours)"

echo "will look for observation at time $OBS_YEAR $OBS_MONTH $OBS_DAY $CPL_SECONDS (seconds)"
echo "will look for observation at time $OBS_YEAR $OBS_MONTH $OBS_DAY $CPL_HOUR (hours)"
#-------------------------------------------------------------------------
# Create temporary working directory for the assimilation and go there
#-------------------------------------------------------------------------

set temp_dir = assimilate_ice
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of ICE.
#-----------------------------------------------------------------------------
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames

if (! -d ${BASEOBSROOT}/) then
   echo "CICE-DART requires daily obs_seq files in a base observation directory."
   echo "The directory ${BASEOBSROOT} is not found.  Exiting"
   exit 2
endif

set OBSFNAME = `printf obs_seq.%04d-%02d-%02d-%05d ${OBS_YEAR} ${OBS_MONTH} ${OBS_DAY} ${CPL_SECONDS}`
set OBS_FILE = ${BASEOBSROOT}/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file ${OBS_FILE}"
   echo "ERROR ... no observation file ${OBS_FILE}"
   exit 2
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit 2
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.
# Since the ensemble manager is not used by dart_to_cice,
# it is OK to set it here and have it used by all routines.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" \
          input.nml.$$ >! input.nml || exit 3
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: Stage the files needed for SAMPLING ERROR CORRECTION and QCEFF
#
# The sampling error correction is a lookup table.
# The tables were originally in the DART distribution, but should
# have been staged to $CASEROOT at setup time.
# Each ensemble size has its own (static) file.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set  MYSTRING = `grep 'sampling_error_correction' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`

if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${CASEROOT}/sampling_error_correction_table.nc
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file found."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit 2
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

# Make this better later but for now (MW2024)
set QCEFF_FILE = ${CASEROOT}/cice_bounded_qceff_table.csv
${COPY} ${QCEFF_FILE} .
#=========================================================================
# Block 3: DART_INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(:) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflation_input',  'posterior_inflation_input',
# inf_out_file_name           = 'prior_inflation_output', 'posterior_inflation_output',
# inf_diag_file_name          = 'prior_obs_infl_out',     'posterior_obs_infl_out',
#
# NOTICE: the archiving scripts require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
#
# The inflation file is essentially a duplicate of the DART model state ...
# For the purpose of this script, they are the output of a previous assimilation,
# so they should be named something like prior_inflate_output.YYYY-MM-DD-SSSSS
#
# NOTICE: inf_initial_from_restart and inf_sd_initial_from_restart are somewhat
# problematic. During the bulk of an experiment, these should be TRUE, since
# we want to read existing inflation files. However, the first assimilation
# might need these to be FALSE and then subsequently be set to TRUE.
# There is now only one way to handle this:
#
# 1) create a cookie file called RUNDIR/cice_inflation_cookie
#    The existence of this file will cause this script to set the
#    namelist appropriately. This script will 'eat' the cookie file
#    to prevent this from happening for subsequent executions. If the
#    inflation file does not exist for them, and it needs to, this script
#    should die. The CESM_DART_config script automatically creates a cookie
#    file to support this option.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#=========================================================================

set  MYSTRING = `grep 'inf_flavor' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep 'inf_initial_from_restart' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
set  POSTE_TF = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e ../cice_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set PRIOR_TF = FALSE

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else

      # Look for inflation files from the previous assimilation
      # This is really ugly- sorry.
      
      # Checking for a prior inflation mean file to use
      
      (ls -rt1 ../cice.output_priorinf_mean.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_mean.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation MEAN file."
         echo "ERROR: expected something like ../cice.output_priorinf_mean.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

      # Checking for a prior inflation sd file to use

      (ls -rt1 ../cice.output_priorinf_sd.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_sd.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation SD file."
         echo "ERROR: expected something like ../cice.input_priorinf_sd.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

   endif
else
   echo "Prior Inflation           not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ( $POSTE_INF > 0 ) then

   if ($POSTE_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."

   else if ( -e ../cice_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set POSTE_TF = FALSE

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else

      # Look for inflation files from the previous assimilation
      # This is really ugly- sorry.

      # Checking for a posterior inflation mean file to use
      
      (ls -rt1 ../cice.output_postinf_mean.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_postinf_mean.nc
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation MEAN file."
         echo "ERROR: expected something like ../cice.output_postinf_mean.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

      # Checking for a posterior inflation sd file to use

      (ls -rt1 ../cice.output_postinf_sd.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_postinf_sd.nc
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation SD file."
         echo "ERROR: expected something like ../cice.output_postinf_sd.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

# Eat the cookie regardless
${REMOVE} ../cice_inflation_cookie

#=========================================================================
# Block 4: Create a set of restart files before DART has modified anything.
#
#   filter has the ability to directly modify the cice restart files
#   i.e. it creates the posterior IN-PLACE.
#   We usually want a prior estimate, so we have to save a copy of the
#   input files before we feed them to filter. If we saved every
#   restart, the directory gets polluted pretty fast, so we overwrite 
#   the same filenames over and over. The timestamps IN the file can 
#   confirm the valid time of the model state.
#
#   At this time we also create a list of files we want to read/modify.
#=========================================================================

echo "`date` -- BEGIN CREATING SAFETY FILES"

# create the list of restart files by dereferencing the pointer files. 
# While we are at it, we have to account for the fact they are 1 dir up.
# CP the input (prior) CICE states for restart purposes.
# The original CICE files will be updated directly by filter.

${REMOVE} cice_restarts.txt

set member = 1
while ( ${member} <= ${ensemble_size} )

   set  SAFETY_FILE = `printf cice_prior.r.%04d.nc ${member}`
   set POINTER_FILE = `printf ../rpointer.ice_%04d ${member}`

   set MYFILE = `head -n 1 $POINTER_FILE`
   set ICE_FILENAME = `echo $MYFILE:t`

   echo "../"${ICE_FILENAME} >> cice_restarts.txt

   ${COPY} ../${ICE_FILENAME} ${SAFETY_FILE} &

   @ member++
end

wait

echo "`date` -- END CREATING SAFETY FILES for all ${ensemble_size} members."

#=========================================================================
# If do parameter estimation, call the followsing block
#=========================================================================
# if ($PARAMETER_ESTIMATION == TRUE ) then
# 
# set member = 1
# while ( $member <= ${ensemble_size} )
#     set inst_string = `printf       _%04d $member`
#     set  member_dir = `printf member_%04d $member`
#  
#     if (! -d ${member_dir}) mkdir ${member_dir}
#         cd ${member_dir}
#  
#         set ICE_FILENAME = `head -n $member ../cice_restarts.txt | tail -n 1`
#         set PARAMETER_FILENAME = "CICE_PARAMETERS_PRIOR"${inst_string}"."${CPL_DATE_EXT}".nc"
#  
#         ln -sf ../${ICE_FILENAME} cice_restart.nc
#  
#         ln -sf ../input.nml .
#  
#         ln -sf ../../${PARAMETER_FILENAME} parameter_prior.nc
#         echo "starting cice_to_dart for member ${member} at "`date`
#         ${EXEROOT}/cice_to_dart >! output.${member}.cice_to_dart  &
#         cd ..
#    @ member ++
#  end
# 
# wait
# 
# endif

set member = 1
while ( $member <= ${ensemble_size} )
   set inst_string = `printf  _%04d $member`
   set  member_dir = `printf member_%04d $member`
   if (! -d ${member_dir}) mkdir ${member_dir}
   cd ${member_dir}
   set ICE_FILENAME = `head -n $member ../cice_restarts.txt | tail -n 1`
   if ( $PARAMETER_ESTIMATION == "TRUE" ) then 
      foreach PARAM ($cice_parameters)
         if ( $CONTINUE_RUN == "FALSE" || $use_mean == "TRUE" ) then   
            ncks -v  fsnow ../$ICE_FILENAME -O temp.nc  
       
            ncrename -v fsnow,$PARAM temp.nc
       
            set par_value = `grep "$PARAM" "${CASEROOT}/user_nl_cice${inst_string}" | cut -d'=' -f2 | sed "s/^ *//"`
            echo ${PARAM} $par_value
            ncap2 -s "${PARAM}[nj,ni]=$par_value" temp.nc temp2.nc
            echo "ncap2 done"
            ncks -v ${PARAM} temp2.nc -A ../$ICE_FILENAME 
            rm -rf temp.nc temp2.nc         
         else   # use filters delta r_snow
            set date_yesterday = `date -d "${CPL_DATE[1]}-${CPL_DATE[2]}-${CPL_DATE[3]} 1 day ago" +%F`-${CPL_DATE[4]}
            set pre_restart = `printf $RUNDIR/${CASE}.cice_%04d.r.${date_yesterday}.nc $member`
            ncks -v $PARAM $pre_restart -A ../$ICE_FILENAME
         endif 
      end
   endif

   ln -sf ../${ICE_FILENAME} cice_restart.nc
   ln -sf ../input.nml .
   echo "input.nml linked"
   @ member ++
   cd ..
end

#=========================================================================
# Block 5: Actually run the assimilation.
#
# >@todo FIXME ... this whole section
#
# REQUIRED DART namelist settings:
# &filter_nml:           async                   = 0,
# &filter_nml:           adv_ens_command         = "no_advance_script",
# &filter_nml:           obs_sequence_in_name    = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name   = 'obs_seq.final'
# &filter_nml:           init_time_days          = -1,
# &filter_nml:           init_time_seconds       = -1,
# &filter_nml:           first_obs_days          = -1,
# &filter_nml:           first_obs_seconds       = -1,
# &filter_nml:           last_obs_days           = -1,
# &filter_nml:           last_obs_seconds        = -1,
#
# &filter_nml: input_restart_file_list  = "cice_restarts.txt"
# &filter_nml: output_restart_file_list = "cice_restarts.txt"
# &filter_nml: output_restarts          = .true.
# &filter_nml: stages_to_write          = 'output'
#=========================================================================

# The cice model_mod.f90:static_init_model() has a hardcoded 'cice.r.nc'
# that must exist. The cice_in,drv_in namelists must also exist in this directory 

set TEMPLATEFILE = `head -n 1 cice_restarts.txt`
ln -sf $TEMPLATEFILE   cice.r.nc
ln -sf ../ice_in_0001  cice_in
ln -sf ../drv_in       drv_in

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit 5
echo "`date` -- END FILTER"

# 1) rename DART files to reflect current date and component
# 2) move to RUNDIR so they get archived and the DART_INFLATION block works next cycle

foreach FILE ( input_*mean.nc   input_*sd.nc   \
               preassim_*nc     \
               postassim_*.nc   \
               output_*mean.nc  output_*sd.nc  \
               dart_log*        obs_seq.final )

   if ( -e $FILE ) then
      set  FEXT = $FILE:e
      set FBASE = $FILE:r
      ${MOVE} $FILE ../cice.${FBASE}.${CPL_DATE_EXT}.${FEXT}
   else
      echo "$FILE does not exist, no need to take action."
   endif

end

# Copy obs_seq.final files to a place that won't be archived,
# so that they don't need to be retrieved from the HPSS.
if (! -d ../../Obs_seqs) mkdir ../../Obs_seqs
${COPY} ../cice.obs_seq.${CPL_DATE_EXT}.final ../../Obs_seqs &

# Handle localization_diagnostics_files
set MYSTRING = `grep 'localization_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set loc_diag = $MYSTRING[2]
if (-f $loc_diag) then
   $MOVE $loc_diag ../cice.${loc_diag}.${CPL_DATE_EXT}
endif

# Handle regression diagnostics
set MYSTRING = `grep 'reg_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set reg_diag = $MYSTRING[2]
if (-f $reg_diag) then
   $MOVE $reg_diag ../cice.${reg_diag}.${CPL_DATE_EXT}
endif

#=========================================================================
# Block 6: 
# The filter settings update the cice netcdf files directly - BUT -
# they need to be rebalanced before being used. The rebalancing is done
# by the dart_to_cice program.
# Each member will do its job in its own directory.
# Block 7: The ice files have now been updated, move them into position.
# >@todo FIXME ... rename 'dart_to_cice' to 'rebalance_cice' or something
# more accurate.
#=========================================================================

echo "`date` -- BEGIN DART-TO-CICE"
set member = 1
while ( $member <= $ensemble_size )

   set inst_string = `printf       _%04d $member`
   set  member_dir = `printf member_%04d $member`

   if (! -d ${member_dir}) mkdir ${member_dir}
   cd ${member_dir}

   ${REMOVE} output.${member}.dart_to_ice

   set ICE_FILENAME = `head -n $member ../cice_restarts.txt | tail -n 1`

   ${LINK} ../${ICE_FILENAME} dart_restart.nc || exit 6

   #========================================================
   # FEI: link the prior restart file cice_prior.r.$member.nc 
   # to cice_restart.nc  `
   #========================================================

   set PRIOR_FILENAME = `printf cice_prior.r.%04d.nc $member`
   ln -sf ../${PRIOR_FILENAME} cice_restart.nc

   # create a postprocessed_restart.nc file for changes to be written to
   cp ../${ICE_FILENAME} postprocessed_restart.nc

   #========================================================
   # FEI: link yesterday's restart file 
   # to pre_restart.nc  
   # MOLLY: link original restart file to original_cice_restart_file
   #========================================================
   # cd ${CASEROOT}

   # if ( ${CONTINUE_RUN} == "TRUE" ) then

   #    set date_yesterday = `date -d "${CPL_DATE[1]}-${CPL_DATE[2]}-${CPL_DATE[3]} 1 day ago" +%F`-${CPL_DATE[4]}
   #    echo "the forecast began at " $date_yesterday
   #    set pre_restart = `printf $RUNDIR/${CASE}.cice_%04d.r.${date_yesterday}.nc $member`

   # else

   #    set infile   = `printf $RUNDIR/ice_in_%04d $member`
   #    set MYSTRING = `grep 'ice_ic' $infile`
   #    set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
   #    set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
   #    set pre_restart = $MYSTRING[2]

   # endif

   # cd $RUNDIR/assimilate_ice/${member_dir}          #go back to the assim dir
#   ${LINK} $pre_restart pre_restart.nc


   echo "starting dart_to_ice for member ${member} at "`date`
   ${EXEROOT}/dart_to_cice >! output.${member}.dart_to_ice 

   # copy the postprocessed_restart.nc back to the restart file in the model run directory
   cp postprocessed_restart.nc ../${ICE_FILENAME}

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.dart_to_ice | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'dart_to_cice' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_cice' ... ERROR"
   exit 2
endif

echo "`date` -- END DART-TO-CICE for all ${ensemble_size} members."


#=========================================================================
# Fei: Now take an average of the 2D parameter field and get a global constant
# Put the new parameter to user_nl_cice_xxxx 
#=========================================================================
if ( $PARAMETER_ESTIMATION == "TRUE" ) then
set member = 1

set npar = `echo $#cice_parameters`

while ( $member <= $ensemble_size )

   set inst_string = `printf       _%04d $member`
   set  member_dir = `printf member_%04d $member`

   if (! -d ${member_dir} ) exit 

   cd ${member_dir}

   set ICE_FILENAME = `head -n $member ../cice_restarts.txt | tail -n 1`

   #${LINK} ../${ICE_FILENAME} dart_restart.nc || exit 6

   #set PARAMETER_FILENAME = "CICE_PARAMETERS_POSTERIOR"${inst_string}"."${CPL_DATE_EXT}".nc"
   set ipar = 1
   while ( $ipar <= $npar )   
      set parameter = $cice_parameters[$ipar]
 
      # use ncwa to get the spatial averaged parameter
      # store in a netcdf file $parameter.nc
     # ncks -v $parameter dart_restart.nc $PARAMETER_FILENAME
   
      cp dart_restart.nc bob.nc 
      ncks -v gw $WORK/headfiles/gw.nc -A bob.nc 
      ncwa -a nj,ni -w gw -v $parameter bob.nc ${parameter}.nc

      rm -rf bob.nc
      # dump the parameter value into a ascii file
      ncdump  -v $parameter ${parameter}.nc > ${parameter}.txt

      # extract the parameter value from the ascii file
      set par_string = `grep -oP "${parameter} =\K.*" ${parameter}.txt`
      set par_value = $par_string[1]

      # now put the new parameter value into user_nl_cice_xxxx under the case root
      set currdir = `pwd`
      cd $CASEROOT
      set myfile = `printf user_nl_cice${inst_string}`
      sed -e "/${parameter}/c\${parameter} = ${par_value} " $myfile >${myfile}.new
      mv ${myfile}.new ${myfile}
      cd ${currdir} 
      @ ipar ++
    end

    cd ..
    @ member ++
 end

endif 
#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END ICE_ASSIMILATE"

exit 0