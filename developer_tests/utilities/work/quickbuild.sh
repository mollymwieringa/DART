#!/usr/bin/env bash

main() {


[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildfunctions.sh

MODEL="template"
LOCATION="threed_sphere"
dev_test=1
TEST="utilities"

serial_programs=(
PrecisionCheck
error_handler_test
file_utils_test
find_enclosing_indices_test
find_first_occurrence_test
nml_test
parse_args_test
sort_test
stacktest
)


# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# preprocess not needed for these tests

# build 
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
