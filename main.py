# This is the main script for applying the SBM to an inlet with large gas bubbles entering a continuous liquid.
# The modelling itself is done in two Python-scripts: the first
# one ("read_inlet_<solver>.py") is flow solver-dependent and is used to create a table from the inlet geometry of the
# case. The second script ("inlet_modeling.py") is solver-independent and creates the actual inlet. 
# The values for the boundary condition are then written by running the third script ("write_bc_<solver>.py")

from utils import os, yaml, shutil, logging, get_openfoam_type
import cProfile
logger = logging.getLogger(__name__)

def main(config):
    # define config variables
    cfd_program = config["packages"]["cfd_program"]
    cfd_version = config["packages"]["cfd_version"]
    config["openfoam_type"] = get_openfoam_type(cfd_version)
    case_path = config["case_path"]
    purge_boundary_data = config["purge_boundary_data"]

    # remove previous log file
    log_path = os.path.join(case_path, "sbm.log")
    if os.path.exists(log_path):
        os.remove(log_path)

    logging.basicConfig(filename=log_path, level=logging.INFO)
    logger.info("++++++++++Synthetic Bubble Model start++++++++++")

    # check if cfd program is defined
    if (cfd_program.lower() != "openfoam") & (cfd_program.lower() != "ansys_cfd"):
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")
    
    if purge_boundary_data:
        boundary_data_path = os.path.join(case_path, "constant", "boundaryData")
        if os.path.exists(boundary_data_path):
            shutil.rmtree(boundary_data_path)

    # =========================================================================================================================================================================== 

    # read SBM parameters
    parameter_path = os.path.join(case_path, "parameter.yaml")
    if not os.path.exists(parameter_path):
        raise RuntimeError(f"parameter.yaml file not found in case directory \n \t {case_path}")
    else:
        with open(parameter_path, "r") as f:
            param = yaml.load(f, Loader=yaml.SafeLoader)

    # Create postProcess function configuration for "writeCellCentres" and "writeCellAreas"
    if cfd_program.lower() == "openfoam":
        from src.create_postprocess_functions import create_postprocess_functions
        create_postprocess_functions(config)

    # Read CFD inlet data 
    if cfd_program.lower() == "openfoam":
        from src.read_inlet_openfoam import read_inlet_openfoam
        read_inlet_openfoam(config, param)
    elif cfd_program.lower() == "ansys_cfd":
        from src.read_inlet_ansyscfd import read_inlet_ansyscfd
        read_inlet_ansyscfd(config, param)
    else:
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    # calculate inlet boundary conditions
    from src.inlet_modelling import inlet_modelling
    inlet_modelling(config, param)

    # Write output
    if cfd_program.lower() == "openfoam":
        from src.write_bc_openfoam import write_bc_openfoam
        write_bc_openfoam(config, param)
    elif cfd_program.lower() == "ansys_cfd":
        from src.write_bc_ansyscfd import write_bc_ansyscfd
        write_bc_ansyscfd(config, param)
    else:
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    logger.info("++++++++++Synthetic Bubble Model end++++++++++")

if __name__=='__main__':

    # read configuration file
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    profile_run = config.get("profile_run", False)

    if profile_run:
        cProfile.run("main(config)", "sbm.prof")
        print(f"Profile results printed to sbm.prof. Visualize by running 'snakeviz sbm.prof'")
    else:
        main(config)