import yaml


if __name__=="__main__":
    print("++++++++++Synthetic Bubble Model start++++++++++")

    # read configuration file
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    cfd_program = config["packages"]["cfd_program"]
    case_path = config["case_path"]

    # check if cfd program is defined
    if (cfd_program.lower() != "openfoam") | (cfd_program.lower() != "ansys_cfd"):
        RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS")


    # =======================================================================================    

    # Read CFD inlet data 
    if cfd_program.lower() == "openfoam":
        from src.read_inlet_openfoam import read_inlet_openfoam
        read_inlet_openfoam(config)
    elif cfd_program.lower() == "ansys_cfd":
        from src.read_inlet_ansyscfd import read_inlet_ansyscfd
        read_inlet_ansyscfd(config)
    else:
        RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS")

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
        RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS")

    print("++++++++++Synthetic Bubble Model end++++++++++")