#!/bin/sh
# This is the main script for applying the SBM to an inlet with large gas bubbles entering a continuous liquid. This
# script is only used to concentrate the user-input for the model (and is dependent on the flow solver the user wants to
# use, denoted by the main IF-statement in the script). In this example Fluent is used. The modelling itself is done in
# two Python-scripts: the first one ("read_inlet_<solver>.py") is flow solver-dependent and is used to create a table
# from the inlet geometry of the case. The second script ("InletModelling.py") is solver-independent and creates the
# actual inlet. The used bubble shapes are defined in "InletModelling.py".

# User input 
export CFD_PROGRAMME=ANSYS_CFD
export CFD_VERSION=2023R1 # Version of the CFD-solver (module name= CFD_PROGRAMME/CFD_VERSION)
export PYTHON_VERSION=Anaconda3-python
export DIM=3 # Number of geometrical dimensions of the case (either 2 or 3)
export CASE_PATH=/cfdfile2/data/fm/henri/GO-VIKING/Box_2 # Location of base case
export startTime=5 # First flow time to be defined
export endTime=10 # Last flow time to be defined
export timeStepSize=0.0001 # Time step size to be used in following calculation
export tunit=0.5 # Unit time scale - parameter for inlet modelling
export inletName=inlet # Name of the inlet boundary to be modelled
export rhog=1.225 # Density of the gas
export rhol=998.2 # Density of the liquid
export mg=0.214375 # Amount of the gas to be introduced in domain over a time 'tunit'
export tol_mg=1e-05 # Tolerance on the amount of gas to be introduced in domain over a time 'tunit'
export U=1 # Velocity of the mixture  to be introduced in domain
export intersectBoundary=True # Boolean indicating whether bubbles may intersect domain boundaries
export intersectBubble=False # Boolean indicating whether new bubble may intersect previously defined bubbles
export mnpf=4 # max nodes per face (Fluent)
export time_step_start=50000 # start_time_step (Fluent)
export n_time_steps=50000 # number of time steps (Fluent)
export case_name=second_test.cas.h5 # name of case file (Fluent)

# Load Python (inlet modelling performed in Python Anaconda)
module load $PYTHON_VERSION

# Execution of the script depends on CFD-programme to be used
export CFD_MODULE=$CFD_PROGRAMME/$CFD_VERSION
# Read the inlet geometry: This script is designed to be case-dependent
if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
then
	python readInlet_"$(echo $CFD_PROGRAMME | tr '[:upper:]' '[:lower:]')".py $DIM $CASE_PATH $case_name $CFD_MODULE $startTime $inletName
elif [ "$CFD_PROGRAMME" = "ANSYS_CFD" ]
then
  python read_inlet_fluent.py $DIM $CASE_PATH $case_name $time_step_start $n_time_steps $inletName $mnpf
else
        echo "The CFD-programme ($CFD_PROGRAMME) you have tried to use is not defined in the main script."
	exit 1
fi

# Model the inlet : This script is designed to be case-independent
python inletModelling.py $CASE_PATH $startTime $endTime $timeStepSize $tunit $inletName $rhog $rhol $mg $tol_mg $U $intersectBoundary $intersectBubble

# The modelled inlet has to be written in a format compatible with the flow solver that is to be used
if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
then
        # Write the inlet boundary condition with Python-script
	python writeBC_"$(echo $CFD_PROGRAMME | tr '[:upper:]' '[:lower:]')".py $CASE_PATH $startTime $inletName
elif [ "$CFD_PROGRAMME" = "ANSYS_CFD" ]
then
  python write_bc_fluent.py $DIM $CASE_PATH $case_name $time_step_start $inletName
else
        echo "The CFD-programme ($CFD_PROGRAMME) you have tried to use is not defined in the main script."
        exit 1
fi
