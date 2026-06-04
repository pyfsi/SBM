# This is the main script for applying the SBM to an inlet with large gas bubbles entering a continuous liquid.
# The modelling itself is done in two Python-scripts: the first
# one ("read_inlet_<solver>.py") is flow solver-dependent and is used to create a table from the inlet geometry of the
# case. The second script ("inlet_modeling.py") is solver-independent and creates the actual inlet. 
# The values for the boundary condition are then written by running the third script ("write_bc_<solver>.py")

from utils import os, yaml, shutil

if __name__=="__main__":
    print("++++++++++Synthetic Bubble Model start++++++++++")

    # read configuration file
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    cfd_program = config["packages"]["cfd_program"]
    case_path = config["case_path"]
    purge_boundary_data = config["purge_boundary_data"]

    # check if cfd program is defined
    if (cfd_program.lower() != "openfoam") & (cfd_program.lower() != "ansys_cfd"):
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")
    
    if purge_boundary_data:
        boundary_data_path = os.path.join(case_path, "constant", "boundaryData")
        if os.path.exists(boundary_data_path):
            shutil.rmtree(boundary_data_path)

    # ======================================================================================= 

    # Create postProcess function configuration for "writeCellCentres" and "writeCellAreas"
    from src.create_postprocess_functions import create_postprocess_functions
    create_postprocess_functions(config)

    # Read CFD inlet data 
    if cfd_program.lower() == "openfoam":
        from src.read_inlet_openfoam import read_inlet_openfoam
        read_inlet_openfoam(config)
    elif cfd_program.lower() == "ansys_cfd":
        from src.read_inlet_ansyscfd import read_inlet_ansyscfd
        read_inlet_ansyscfd(config)
    else:
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    # calculate inlet boundary conditions
    from src.inlet_modelling import inletModelling
    inletModelling(config)

    # Write output
    if cfd_program.lower() == "openfoam":
        from src.write_bc_openfoam import write_bc_openfoam
        write_bc_openfoam(config)
    elif cfd_program.lower() == "ansys_cfd":
        from src.write_bc_ansyscfd import write_bc_ansyscfd
        write_bc_ansyscfd(config)
    else:
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    print("++++++++++Synthetic Bubble Model end++++++++++")