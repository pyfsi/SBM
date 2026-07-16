#!/bin/sh
test=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# define path to SBM directory
echo "Defining paths SBM and SBM_POSTPROC"
export SBM=$test
export SBM_POSTPROC=$test/postprocess

# append SBM directory path to PYTHONPATH
echo "Adding $test to PYTHONPATH"
export PYTHONPATH=$test:$PYTHONPATH
