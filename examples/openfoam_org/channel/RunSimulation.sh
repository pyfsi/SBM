#!/bin/bash

echo "Creating mesh for geometry"
blockMesh
echo ""

echo "Running SBM script"
python $SBM/main.py
echo ""

echo "Starting OpenFOAM simulation"
bash ./Allrun
echo ""

echo "RunSimulation.sh ended"
