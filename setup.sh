#!/bin/sh
sbm=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# define path to SBM directory
echo "Defining paths SBM and SBM_POSTPROC"
export SBM=$sbm
export SBM_POSTPROC=$sbm/postprocess

# append SBM directory path to PYTHONPATH
echo "Adding $sbm to PYTHONPATH"
export PYTHONPATH=$sbm:$PYTHONPATH
