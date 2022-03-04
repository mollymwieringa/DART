#!/usr/bin/env bash

main() {

[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildconvfunctions.sh

CONVERTER=ok_mesonet
LOCATION=threed_sphere
EXTRA="$DART/observations/obs_converters/obs_error/dewpoint_obs_err_mod.f90
       $DART/observations/obs_converters/obs_error/ncep_obs_err_mod.f90
       $DART/observations/obs_converters/MADIS/meteor_mod.f90"


programs=(
convert_ok_mesonet
read_geo
obs_sequence_tool
advance_time
)

# build arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildconv


# clean up
\rm -f -- *.o *.mod

}

main "$@"
