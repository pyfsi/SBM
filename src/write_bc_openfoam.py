from utils import np, sys, os, logging, shutil, subprocess, ThreadPoolExecutor
from utils import truncate
from utils import SBM_OUTPUT
from functools import partial
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

def write_time_step(i, args):
    # arguments
    delta_time = args[0]
    boundary_inlet_path = args[1]
    alpha_name = args[2]
    n_inlet_points = args[3]
    U_inlet = args[4]
    VOF_inlet = args[5]
    time_inlet = args[6]

    time_i = time_inlet[i]
    trunc_time_i = truncate(time_i, delta_time)

    # Set location of files
    U_out_path = f"{boundary_inlet_path}/{trunc_time_i}/U"
    VOFw_out_path = f"{boundary_inlet_path}/{trunc_time_i}/{alpha_name}"

    # make directory for each time step
    time_i_dir = os.path.join(boundary_inlet_path, str(trunc_time_i))
    os.mkdir(time_i_dir)
    make_u_dir = f"touch {boundary_inlet_path}/{trunc_time_i}/U;"
    subprocess.run(make_u_dir, shell=True)
    make_alpha_dir = f"touch {boundary_inlet_path}/{trunc_time_i}/{alpha_name};"
    subprocess.run(make_alpha_dir, shell=True)

    # Write OpenFOAM-header
    write_header(U_out_path, "vectorAverageField", "values")
    write_header(VOFw_out_path, "scalarAverageField", "values")

    # write content
    with open(U_out_path, 'a+') as f:
        f.write('//Data points'+'\n')
        f.write(str(n_inlet_points)+'\n'+'('+'\n')
        for j in np.arange(n_inlet_points):
            f.write('(' + str(U_inlet[j, i, 0]) + ' ' + str(U_inlet[j, i, 1]) + ' ' + str(U_inlet[j, i, 2]) + ') \n')
        f.write(')'+'\n')

    with open(VOFw_out_path, 'a+') as f:
        f.write('//Data points'+'\n')
        f.write(str(n_inlet_points)+'\n'+'('+'\n')
        for j in np.arange(n_inlet_points):
            f.write(str(VOF_inlet[j,i,0])+'\n')
        f.write(')'+'\n')

    # Write OpenFOAM-footer
    write_footer(U_out_path)
    write_footer(VOFw_out_path)

def write_bc_openfoam(config):
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
    command = f"grep -nrw {inlet_name} {start_time}/U | cut -d : -f 1 > lineNr_U"
    subprocess.run(command, cwd=case_path, shell=True)
    lineNrU_path = os.path.join(case_path, "lineNr_U")
    lineNameNr = int(open(lineNrU_path, 'r').readline())
    lineTypeU = lineNameNr+1  # inlet BC type is defined on this line
    U_path = os.path.join(case_path, str(start_time), "U")
    readTypeU = (open(U_path).readlines())[lineTypeU]
    boundaryConditionU = (readTypeU.split())[-1][0:-1]
    os.remove(lineNrU_path)
    if boundaryConditionU != "timeVaryingMappedFixedValue":
        logger.error(f"U at {inlet_name} has type {readTypeU}")
        sys.exit(f"The condition for 'U' at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \n")

    # Boundary condition should not be setAverage
    command = f"grep -nr 'setAverage' {start_time}/U | cut -d : -f 1 > lineNr_setAvg"
    subprocess.run(command, cwd=case_path, shell=True)
    lineNrsetAvg_path = os.path.join(case_path, "lineNr_setAvg")
    lineNameNr = int(open(lineNrsetAvg_path, 'r').readline())-1  # Line where 'setAverage' is defined
    readSetAvg = (open(U_path).readlines())[lineNameNr]
    setAvg = (readSetAvg.split())[-1][0:-1]
    os.remove(lineNrsetAvg_path)
    if setAvg != "false":
        logger.error(f"Boundary condition at {inlet_name} should not contain setAverage for 'U'")
        sys.exit(f"The boundary condition at '{inlet_name}' defines an averaging operation for 'U'. \
                This is not compatible with the transient inlet modelling defined in the Python script. \n")

    # Check boundary condition for alpha
    command = f"grep -nrw {inlet_name} {start_time}/{alpha_name} | cut -d : -f 1 > lineNr_VOFw"
    subprocess.run(command, cwd=case_path, shell=True)
    lineNr_VOFw_path = os.path.join(case_path, "lineNr_VOFw")
    lineNameNr = int(open(lineNr_VOFw_path, 'r').readline())
    lineTypeVOFw = lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
    alpha_path = os.path.join(case_path, str(start_time), alpha_name)
    readTypeVOFw = (open(alpha_path).readlines())[lineTypeVOFw]
    boundaryConditionVOFw = (readTypeVOFw.split())[-1][0:-1]
    os.remove(lineNr_VOFw_path)
    if boundaryConditionVOFw != "timeVaryingMappedFixedValue":
        sys.exit(f"The condition for {alpha_name} at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \
                See line {lineNameNr} \n")
    command = f"grep -nr 'setAverage' {start_time}/{alpha_name} | cut -d : -f 1 > lineNr_setAvg"
    subprocess.run(command, cwd=case_path, shell=True)
    lineNr_setAvg_path = os.path.join(case_path, "lineNr_setAvg")

    try:
        lineNameNr = int(lineNr_setAvg_path)-1 #Line where 'setAverage' is defined
        readSetAvg = (open(alpha_path).readlines())[lineNameNr]
        setAvg = (readSetAvg.split())[-1][0:-1]
        os.remove(lineNr_setAvg_path)
        if setAvg != "false":
            sys.exit(f"Error! The boundary condition at '{inlet_name}' defines an averaging operation for {alpha_name}. \
                    This is not compatible with the transient inlet modelling defined in the Python script.")
        print(f"Inlet definition in folder {start_time} is OK.")
    except:
        os.remove(lineNr_setAvg_path)

    # Prepare OpenFOAM-directory
    logger.info("Creating folder 'boundaryData'.")
    constant_path = os.path.join(case_path, "constant")
    if not os.path.exists(constant_path):
        raise RuntimeError("'constant' folder does not exist in the OpenFOAM case directory.")
    boundary_data_path = os.path.join(case_path, "constant", "boundaryData")
    if os.path.exists(boundary_data_path):
        raise RuntimeError("'boundaryData' folder already exists.")

    # Create boundaryData
    os.mkdir(os.path.join(case_path, "constant", "boundaryData"))
    os.mkdir(os.path.join(case_path, "constant", "boundaryData", inlet_name))
    boundary_inlet_path = os.path.join(case_path, "constant", "boundaryData", inlet_name)
    command = f"touch {boundary_inlet_path}/points"
    subprocess.run(command, shell=True)

    # Write 'points'-file
    points_out_path = f"{case_path}/constant/boundaryData/{inlet_name}/points"
    n_inlet_points = len(inlet_coord[:, 0])
    with open(points_out_path, 'a+') as f:
        f.write(str(n_inlet_points)+'\n')
        f.write('('+'\n')
        for j in np.arange(n_inlet_points):
            f.write('(' + str(inlet_coord[j, 1]) + ' ' + str(inlet_coord[j, 2]) + ' ' + str(inlet_coord[j, 3]) + ') \n')
        f.write(')'+'\n')

    # Write velocity and volume fraction for water at each time step.
    logger.info("Writing boundary condition to folder 'boundaryData'.")
    args = [delta_time, boundary_inlet_path, alpha_name, n_inlet_points, U_inlet, VOF_inlet, time_inlet]
    partial_write = partial(write_time_step, args=args)
    n_workers = min(8, os.cpu_count())
    logger.info(f"Using {n_workers} threads for writing routine.")
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        executor.map(partial_write, np.arange(len(time_inlet)))
    logger.info("Boundary condition was successfully saved in 'boundaryData'.")

    # move area postProcess file
    target_area_path = os.path.join(output_path, "area")
    if os.path.exists(target_area_path):
        os.remove(target_area_path)
    area_path = os.path.join(case_path, "0", "area")
    if os.path.exists(area_path):
        shutil.move(area_path, output_path)

    logger.info("Finished writing boundary condition to folder 'boundaryData'.")
