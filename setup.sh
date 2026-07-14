#!/bin/sh

# get path to sbm directory
export SBM=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# load the specified Python module
export PYTHON_VERSION=Anaconda3-python/2023.09-0
module load $PYTHON_VERSION
