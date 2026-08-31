#!/bin/bash

echo "Creating mesh for geometry"
cd "../meshCase"
bash ./Allmesh
cd "../case"

echo "Copying mesh from meshCase and renumbering patches"
bash ./Allpreprocess

echo "Running SBM script"
python $SBM/main.py

echo "Starting OpenFOAM simulation"
bash ./Allrun
