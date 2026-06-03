#!/bin/sh

# load the specified Python module
export PYTHON_VERSION=Anaconda3-python/2023.09-0
module load $PYTHON_VERSION

# Run Python without creating pycache fils
PYTHONDONTWRITEBYTECODE=1 python main.py -p no:cacheprovider
