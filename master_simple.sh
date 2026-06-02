#!/bin/sh
# This is the main script for applying the SBM to an inlet with large gas bubbles entering a continuous liquid. This
# script is only used to concentrate the user-input for the model (and is dependent on the flow solver the user wants to
# use, denoted by the main IF-statement in the script). The modelling itself is done in two Python-scripts: the first
# one ("read_inlet_<solver>.py") is flow solver-dependent and is used to create a table from the inlet geometry of the
# case. The second script ("InletModelling.py") is solver-independent and creates the actual inlet. The used bubble
# shapes are defined in "InletModelling.py".

# 
export CFD_PROGRAMME=OpenFOAM # OpenFOAM or ANSYS_CFD
export CFD_VERSION=v2312-foss-2023a # OpenFOAM/v2406-foss-2023a Version of the CFD-solver (module name= CFD_PROGRAMME/CFD_VERSION)
export PYTHON_VERSION=Anaconda3-python/2023.09-0 #Anaconda3-python/2022.10
export CASE_PATH=/cfdfile1/data/fm/radiputr/Documents/02_Simulation/SBM/examples/channel_bubble # Location of base case

# Load Python (inlet modelling performed in Python Anaconda)
module load $PYTHON_VERSION

# Execution of the script depends on CFD-programme to be used
export CFD_MODULE=$CFD_PROGRAMME/$CFD_VERSION

# remove old boundary data TODO delete this
export OLD_BOUNDARY_DATA=$CASE_PATH/constant/boundaryData
rm -rf $OLD_BOUNDARY_DATA

python main.py

# # Read the inlet geometry: This script is designed to be case-dependent
# if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
# then
# 	# python readInlet_OpenFOAM.py $DIM $CASE_PATH $CFD_MODULE $startTime $inletName
#   python readInlet_OpenFOAM.py
# elif [ "$CFD_PROGRAMME" = "ANSYS_CFD" ]
# then
#   python read_inlet_fluent.py $DIM $CASE_PATH $case_name $time_step_start $n_time_steps $inletName $mnpf $CORES
# else
#   echo "The CFD-programme ($CFD_PROGRAMME) you have tried to use is not defined in the main script."
# 	exit 1
# fi

# # Model the inlet : This script is designed to be case-independent
# # python inletModelling.py $CASE_PATH $startTime $endTime $timeStepSize $tunit $inletName $rhog $rhol $mg $tol_mg $U $intersectBoundary $intersectBubble $mgb_min $mgb_max
# python inletModelling.py

# # The modelled inlet has to be written in a format compatible with the flow solver that is to be used
# if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
# then
# # Write the inlet boundary condition with Python-script
# 	# python writeBC_OpenFOAM.py $CASE_PATH $startTime $inletName $fluid_name
#   python writeBC_OpenFOAM.py
# elif [ "$CFD_PROGRAMME" = "ANSYS_CFD" ]
# then
#   python write_bc_fluent.py $DIM $CASE_PATH $case_name $time_step_start $inletName $CORES
# else
#   echo "The CFD-programme ($CFD_PROGRAMME) you have tried to use is not defined in the main script."
#   exit 1
# fi
