from utils import np, sys, os, logging, shutil, truncate, SBM_OUTPUT
logger = logging.getLogger(__name__)

def write_header(file_loc: str, class_name: str, object_name: str):
    with open(file_loc, 'w') as f:
        f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
        f.write(r'| =========                 |                                                 |'+"\n")
        f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
        f.write(r'|  \\    /   O peration     | Version:  x                                   |'+"\n")
        f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
        f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
        f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
        f.write(r'FoamFile'+"\n")
        f.write(r'{'+"\n")
        # f.write('\t version \t\t x;'+"\n")
        f.write('\t format \t\t ascii;'+"\n")
        f.write('\t class \t\t ' + class_name + ';'+"\n")
        f.write('\t object \t\t ' + object_name + ';'+"\n");
        f.write('}'+"\n")
        f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
        f.write("\n")

def write_footer(file_loc: str):
    with open(file_loc, 'a+') as f:
        f.write("\n")
        f.write(r'// ************************************************************************* //' + "\n")

def write_bc_openfoam(config):
    logger.info("========================Start write_bc_openfoam========================")
    print(f"Running write_bc_openfoam")

    # config variables
    case_path = os.getcwd()
    start_time = int(config["cfd"]["start_time"])
    inlet_name = str(config["cfd"]["inlet_name"])
    alpha_name = "alpha."+config["cfd"]["fluid_name"]
    delta_time = float(config["cfd"]["delta_time"])
    output_path = os.path.join(case_path, SBM_OUTPUT)

    logger.info("Reading inlet values from NPY-files")
    U_inlet = np.load(os.path.join(output_path, "inletDefinition-U.npy"))
    VOF_inlet = np.load(os.path.join(output_path, "inletDefinition-VOFw.npy"))
    time_inlet = np.load(os.path.join(output_path, "inletDefinition-time.npy"))
    inlet_coord = np.load(os.path.join(output_path, "inletPython.npy"))
    logger.info("Finished reading NPY-files.")

    logger.info(f"Checking inlet definition in folder {start_time}.")
    # Check boundary condition for 'U'. Boundary condition should be type timeVaryingMappedFixedValue
    os.system(f"cd {case_path}; grep -nrw {inlet_name} {start_time}/U | cut -d : -f 1 > lineNr_U")
    lineNrU_path = os.path.join(case_path, "lineNr_U")
    lineNameNr = int(open(lineNrU_path, 'r').readline())
    lineTypeU = lineNameNr+1  # inlet BC type is defined on this line
    U_path = os.path.join(case_path, str(start_time), "U")
    readTypeU = (open(U_path).readlines())[lineTypeU]
    boundaryConditionU = (readTypeU.split())[-1][0:-1]
    os.system(f"rm {case_path}/lineNr_U")
    if boundaryConditionU != "timeVaryingMappedFixedValue":
        logger.error(f"U at {inlet_name} has type {readTypeU}")
        sys.exit(f"The condition for 'U' at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \n")

    # Boundary condition should not be setAverage
    os.system(f"cd {case_path}; grep -nr 'setAverage' {start_time}/U | cut -d : -f 1 > lineNr_setAvg" )
    lineNrsetAvg_path = os.path.join(case_path, "lineNr_setAvg")
    lineNameNr = int(open(lineNrsetAvg_path, 'r').readline())-1  # Line where 'setAverage' is defined
    readSetAvg = (open(U_path).readlines())[lineNameNr]
    setAvg = (readSetAvg.split())[-1][0:-1]
    os.system(f"rm {case_path}/lineNr_setAvg")
    if setAvg != "false":
        logger.error(f"Boundary condition at {inlet_name} should not contain setAverage for 'U'")
        sys.exit(f"The boundary condition at '{inlet_name}' defines an averaging operation for 'U'. \
                This is not compatible with the transient inlet modelling defined in the Python script. \n")

    # Check boundary condition for alpha
    os.system(f"cd {case_path}; grep -nrw {inlet_name} {start_time}/{alpha_name} | cut -d : -f 1 > lineNr_VOFw" )
    lineNr_VOFw_path = os.path.join(case_path, "lineNr_VOFw")
    lineNameNr = int(open(lineNr_VOFw_path, 'r').readline())
    lineTypeVOFw = lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
    alpha_path = os.path.join(case_path, str(start_time), alpha_name)
    readTypeVOFw = (open(alpha_path).readlines())[lineTypeVOFw]
    boundaryConditionVOFw = (readTypeVOFw.split())[-1][0:-1]
    os.system(f"rm {case_path}/lineNr_VOFw")
    if boundaryConditionVOFw != "timeVaryingMappedFixedValue":
        sys.exit(f"The condition for {alpha_name} at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \
                See line {lineNameNr} \n")
    os.system(f"cd {case_path}; grep -nr 'setAverage' {start_time}/{alpha_name} | cut -d : -f 1 > lineNr_setAvg" )
    lineNr_setAvg_path = os.path.join(case_path, "lineNr_setAvg")

    try:
        lineNameNr = int(lineNr_setAvg_path)-1 #Line where 'setAverage' is defined
        readSetAvg = (open(alpha_path).readlines())[lineNameNr]
        setAvg = (readSetAvg.split())[-1][0:-1]
        os.system(f"rm {case_path}/lineNr_setAvg")
        if setAvg != "false":
            sys.exit(f"Error! The boundary condition at '{inlet_name}' defines an averaging operation for {alpha_name}. \
                    This is not compatible with the transient inlet modelling defined in the Python script.")
        print(f"Inlet definition in folder {start_time} is OK.")
    except:
        os.system(f"rm {case_path}/lineNr_setAvg")

    # Prepare OpenFOAM-directory
    logger.info("Creating folder 'boundaryData'.")
    constant_path = os.path.join(case_path, "constant")
    if not os.path.exists(constant_path):
        raise RuntimeError("'constant' folder does not exist in the OpenFOAM case directory.")
    boundary_data_path = os.path.join(case_path, "constant", "boundaryData")
    if os.path.exists(boundary_data_path):
        raise RuntimeError("'boundaryData' folder already exists.")

    # Try creating time step directories in 'constant' if it exists
    n_errors = 0
    try:
        n_errors += os.system(f"cd {case_path}/constant; mkdir boundaryData; mkdir boundaryData/{inlet_name}")
        n_errors += os.system(f"touch {case_path}/constant/boundaryData/{inlet_name}/points")

        for time_i in time_inlet:
            trunc_time_i = truncate(time_i, delta_time)

            make_time_dir = f"cd {case_path}/constant/boundaryData/{inlet_name}; mkdir {trunc_time_i};"
            n_errors += os.system(make_time_dir)

            make_u_dir = f"touch {case_path}/constant/boundaryData/{inlet_name}/{trunc_time_i}/U;"
            n_errors += os.system(make_u_dir)

            make_alpha_dir = f"touch {case_path}/constant/boundaryData/{inlet_name}/{trunc_time_i}/{alpha_name};"
            n_errors += os.system(make_alpha_dir)
        if n_errors > 0:
            sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
    except:
        sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
    logger.info("Boundary condition was successfully saved in 'boundaryData'.")

    # Write velocity and alpha water for every time step.
    # Also create 'inlet/points' containing to store position, where velocity and alpha water are stored
    logger.info("Writing inlet definition to folder 'boundaryData'.")
    n_inlet_points = len(inlet_coord[:, 0])
    for i, time_i in enumerate(time_inlet):
        trunc_time_i = truncate(time_i, delta_time)

        # Set location of files
        points_out_path = f"{case_path}/constant/boundaryData/{inlet_name}/points"
        U_out_path = f"{case_path}/constant/boundaryData/{inlet_name}/{trunc_time_i}/U"
        VOFw_out_path = f"{case_path}/constant/boundaryData/{inlet_name}/{trunc_time_i}/{alpha_name}"

        # Write OpenFOAM-header
        write_header(points_out_path, "vectorField", "points")
        write_header(U_out_path, "vectorAverageField", "values")
        write_header(VOFw_out_path, "scalarAverageField", "values")

        # Write 'points'-file
        with open(points_out_path, 'a+') as f:
            f.write(str(n_inlet_points)+'\n')
            f.write('('+'\n')
            for j in np.arange(n_inlet_points):
                f.write('(' + str(inlet_coord[j, 1]) + ' ' + str(inlet_coord[j, 2]) + ' ' + str(inlet_coord[j, 3]) + ') \n')
            f.write(')'+'\n')

        with open(U_out_path, 'a+') as f:
            # necessary for vectorAverageField, but not used if 'setAverage'
            # is set to false (checked earlier in this script) TODO
            # f.write('//Average'+'\n')
            # f.write('(0 0 0)'+'\n'+'\n')

            f.write('//Data points'+'\n')
            f.write(str(n_inlet_points)+'\n'+'('+'\n')
            for j in np.arange(n_inlet_points):
                f.write('(' + str(U_inlet[j, i, 0]) + ' ' + str(U_inlet[j, i, 1]) + ' ' + str(U_inlet[j, i, 2]) + ') \n')
            f.write(')'+'\n')
            f.close()
            # Write alpha-file
            f = open(VOFw_out_path, 'a+')

            # necessary for vectorAverageField, but not used if 'setAverage'
            # is set to false (checked earlier in this script) TODO
            # f.write('//Average'+'\n')
            # f.write('0'+'\n'+'\n')  #

            f.write('//Data points'+'\n')
            f.write(str(n_inlet_points)+'\n'+'('+'\n')
            for j in np.arange(n_inlet_points):
                f.write(str(VOF_inlet[j,i,0])+'\n')
            f.write(')'+'\n')

        # Write OpenFOAM-footer
        write_footer(points_out_path)
        write_footer(U_out_path)
        write_footer(VOFw_out_path)

    # move area postProcess file
    target_area_path = os.path.join(output_path, "area")
    if os.path.exists(target_area_path):
        os.remove(target_area_path)
    area_path = os.path.join(case_path, "0", "area")
    if os.path.exists(area_path):
        shutil.move(area_path, output_path)

    logger.info("Finished writing boundary condition to folder 'boundaryData'.")
    logger.info("========================End write_bc_openfoam========================")
