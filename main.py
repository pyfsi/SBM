# This is the main script for applying the SBM to an inlet with large gas bubbles entering a continuous liquid.
# The modelling itself is done in two Python-scripts: the first
# one ("read_inlet_<solver>.py") is flow solver-dependent and is used to create a table from the inlet geometry of the
# case. The second script ("inlet_modeling.py") is solver-independent and creates the actual inlet.
# The values for the boundary condition are then written by running the third script ("write_bc_<solver>.py")

from utils import os, yaml, shutil, logging
from utils import get_openfoam_type
from utils import SBM_OUTPUT
import cProfile
logger = logging.getLogger(__name__)


def main(config):
    # define config variables
    cfd_program = config["packages"]["cfd_program"]
    cfd_version = config["packages"]["cfd_version"]
    config["openfoam_type"] = get_openfoam_type(cfd_version)
    case_path = os.getcwd()
    purge_boundary_data = config.get("purge_boundary_data", True)

    # remove previous sbm output files
    sbm_out_path = os.path.join(case_path, SBM_OUTPUT)
    if os.path.exists(sbm_out_path):
        shutil.rmtree(sbm_out_path)
        print(f"Deleting previous SBM output files.")
    os.mkdir(sbm_out_path)

    log_path = os.path.join(sbm_out_path, "sbm.log")
    logging.basicConfig(filename=log_path, level=logging.INFO)
    logger.info("++++++++++Synthetic Bubble Model start++++++++++")

    # check if cfd program is defined
    if (cfd_program.lower() != "openfoam") & (cfd_program.lower() != "ansys_cfd"):
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    if purge_boundary_data:
        boundary_data_path = os.path.join(case_path, "constant", "boundaryData")
        if os.path.exists(boundary_data_path):
            shutil.rmtree(boundary_data_path)
    else:
        raise RuntimeError("ERROR: boundaryData directory already exists.")

    # Create postProcess function configuration for "writeCellCentres" and "writeCellAreas"
    if cfd_program.lower() == "openfoam":
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
    from src.inlet_modelling import inlet_modelling
    inlet_modelling(config)

    # Write output
    if cfd_program.lower() == "openfoam":
        from src.write_bc_openfoam import write_bc_openfoam
        write_bc_openfoam(config)
    elif cfd_program.lower() == "ansys_cfd":
        from src.write_bc_ansyscfd import write_bc_ansyscfd
        write_bc_ansyscfd(config)
    else:
        raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    logger.info("++++++++++Synthetic Bubble Model end++++++++++")

if __name__=='__main__':

    # read configuration file
    config_path = os.path.join(os.getcwd(), "config.yaml")
    with open(config_path, "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    profile_run = config.get("profile_run", False)

    if profile_run:
        prof_name = "sbm.prof"
        prof_path = os.path.join(os.getcwd(), SBM_OUTPUT, prof_name)
        cProfile.run("main(config)", prof_path)
    else:
        main(config)

    print(f"SBM script ended succesfully. See 'sbm.log' for more details.")