#!/bin/bash

echo "Creating mesh for geometry"
blockMesh

echo "Running SBM script"
python $SBM/main.py

echo "Starting OpenFOAM simulation"
bash ./Allrun
