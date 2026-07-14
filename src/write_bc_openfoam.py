# This script is developed to read in the output from the flow solver-independent "inlet_modelling" script, search for a
# specific inlet boundary and rewrite the cell centres and faces on that boundary to an array to be used for the inlet
# modelling in a later stage (not in this script). This output contains two matrices, one for flow speed and one for
# Volume-Of-Fluid of water, where the rows are the cell IDs and the columns are the time steps for which U/VOFw are
# defined. These matrices are converted into text files, to be written in the OpenFOAM-directory such that the
# OpenFOAM-utility "timeVaryingMappedFixedValue" can be used. This script is called automatically by the masterscript
# 'TubeBundle_master.sh', so the user input is channeled to this python script from the bash-script directly.

from utils import np, sys, os, logging, shutil, truncate, NPY_OUT_FOLDER
logger = logging.getLogger(__name__)

#writeHeader: to write OpenFOAM-header in file at location 'fileLoc' - class and object of parameter should be given to function
def writeHeader(fileLoc,className,objectName):
    f = open(fileLoc, 'w')
    f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
    f.write(r'| =========                 |                                                 |'+"\n")
    f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
    f.write(r'|  \\    /   O peration     | Version:  4.x                                   |'+"\n")
    f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
    f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
    f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
    f.write(r'FoamFile'+"\n")
    f.write(r'{'+"\n")
    f.write('\t version \t\t 4.1;'+"\n")
    f.write('\t format \t\t ascii;'+"\n")
    f.write('\t class \t\t ' + className + ';'+"\n")
    f.write('\t object \t\t ' + objectName + ';'+"\n");
    f.write('}'+"\n")
    f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
    f.write("\n")
    f.close()


#writeFooter: to write OpenFOAM-footer in file at location 'fileLoc'
def writeFooter(fileLoc):
    f = open(fileLoc, 'a+')
    f.write("\n")
    f.write(r'// ************************************************************************* //' + "\n")
    f.close()


def write_bc_openfoam(config):
    logger.info("========================Start write_bc_openfoam========================")
    print(f"Running write_bc_openfoam")

    casePath = os.getcwd()
    startTime = int(config["cfd"]["start_time"])
    inletName = str(config["cfd"]["inlet_name"])
    alpha_name = "alpha."+config["cfd"]["fluid_name"]
    delta_time = float(config["cfd"]["delta_time"])
    output_path = os.path.join(casePath, NPY_OUT_FOLDER)

    # Inlet model has been previously stored in "{output_path}/inletDefinition-U.npy" and "{output_path}/inletDefinition-VOFw.npy"
    # Time instants at which these variables were defined are stored in "{output_path}/inletDefinition-time.npy"
    # Coordinates of the inlet geometry points are stored in "{output_path}/inletPython.npy"
    logger.info("Reading defined inlet from Python-files")
    # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of velocity
    UVal = np.load(os.path.join(output_path, "inletDefinition-U.npy"))
    # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of VOFw
    VOFwVal = np.load(os.path.join(output_path, "inletDefinition-VOFw.npy"))
    # List containing the time instants where U and alpha are defined
    timeVal = np.load(os.path.join(output_path, "inletDefinition-time.npy"))
    coordList = np.load(os.path.join(output_path, "inletPython.npy"))
    logger.info("Finished reading NPY-files.")

    # This code writes the inlet values assuming the boundary condition 'timeVaryingMappedFixedValue' is used. It is first checked whether this BC is indeed used. Otherwise, the script returns an error.
    # Also, it is checked whether an averaging operation ('setAverage    true;' in OF) is defined - if so, another error is returned.
    logger.info(f"Checking inlet definition in folder {startTime}.")
    # Check boundary condition for 'U'. Boundary condition should be type timeVaryingMappedFixedValue
    os.system(f"cd {casePath}; grep -nrw {inletName} {startTime}/U | cut -d : -f 1 > lineNr_U")
    lineNrU_path = os.path.join(casePath, "lineNr_U")
    lineNameNr = int(open(lineNrU_path, 'r').readline())
    lineTypeU = lineNameNr+1  # inlet BC type is defined on this line - considering Python starts at index zero
    U_path = os.path.join(casePath, str(startTime), "U")
    readTypeU = (open(U_path).readlines())[lineTypeU]
    boundaryConditionU = (readTypeU.split())[-1][0:-1]
    os.system(f"rm {casePath}/lineNr_U")
    if boundaryConditionU != "timeVaryingMappedFixedValue":
        logger.error(f"U at {inletName} has type {readTypeU}")
        sys.exit(f"The condition for 'U' at the boundary '{inletName}' is not set to 'timeVaryingMappedFixedValue'. \n")

    # Boundary condition should not be setAverage
    os.system(f"cd {casePath}; grep -nr 'setAverage' {startTime}/U | cut -d : -f 1 > lineNr_setAvg" )
    lineNrsetAvg_path = os.path.join(casePath, "lineNr_setAvg")
    lineNameNr = int(open(lineNrsetAvg_path, 'r').readline())-1  # Line where 'setAverage' is defined
    readSetAvg = (open(U_path).readlines())[lineNameNr]
    setAvg = (readSetAvg.split())[-1][0:-1]
    os.system(f"rm {casePath}/lineNr_setAvg")
    if setAvg != "false":
        logger.error(f"Boundary condition at {inletName} should not contain setAverage for 'U'")
        sys.exit(f"The boundary condition at '{inletName}' defines an averaging operation for 'U'. \
                This is not compatible with the transient inlet modelling defined in the Python script. \n")

    # Check boundary condition for alpha
    os.system(f"cd {casePath}; grep -nrw {inletName} {startTime}/{alpha_name} | cut -d : -f 1 > lineNr_VOFw" )
    lineNr_VOFw_path = os.path.join(casePath, "lineNr_VOFw")
    lineNameNr = int(open(lineNr_VOFw_path, 'r').readline())
    lineTypeVOFw = lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
    alpha_path = os.path.join(casePath, str(startTime), alpha_name)
    readTypeVOFw = (open(alpha_path).readlines())[lineTypeVOFw]
    boundaryConditionVOFw = (readTypeVOFw.split())[-1][0:-1]
    os.system(f"rm {casePath}/lineNr_VOFw")
    if boundaryConditionVOFw != "timeVaryingMappedFixedValue":
        sys.exit(f"The condition for {alpha_name} at the boundary '{inletName}' is not set to 'timeVaryingMappedFixedValue'. \
                See line {lineNameNr} \n")
    os.system(f"cd {casePath}; grep -nr 'setAverage' {startTime}/{alpha_name} | cut -d : -f 1 > lineNr_setAvg" )
    lineNr_setAvg_path = os.path.join(casePath, "lineNr_setAvg")

    try:
        lineNameNr = int(lineNr_setAvg_path)-1 #Line where 'setAverage' is defined
        readSetAvg = (open(alpha_path).readlines())[lineNameNr]
        setAvg = (readSetAvg.split())[-1][0:-1]
        os.system(f"rm {casePath}/lineNr_setAvg")
        if setAvg != "false":
            sys.exit(f"Error! The boundary condition at '{inletName}' defines an averaging operation for {alpha_name}. \
                    This is not compatible with the transient inlet modelling defined in the Python script.")
        print(f"Inlet definition in folder {startTime} is OK.")
    except:
        os.system(f"rm {casePath}/lineNr_setAvg")

    # Prepare OpenFOAM-directory
    logger.info("Creating folder 'boundaryData'.")
    constant_path = os.path.join(casePath, "constant")
    if not os.path.exists(constant_path):
        raise RuntimeError("'constant' folder does not exist in the OpenFOAM case directory.")
    boundary_data_path = os.path.join(casePath, "constant", "boundaryData")
    if os.path.exists(boundary_data_path):
        raise RuntimeError("'boundaryData' folder already exists.")

    errVal = 0  # Integer denoting whether os.system has error (>0: at least one error)
    try:
        errVal += os.system(f"cd {casePath}/constant; mkdir boundaryData; mkdir boundaryData/{inletName}")
        errVal += os.system(f"touch {casePath}/constant/boundaryData/{inletName}/points")

        for time_i in timeVal:
            trunc_time_i = truncate(time_i, delta_time)

            make_time_dir = f"cd {casePath}/constant/boundaryData/{inletName}; mkdir {trunc_time_i};"
            errVal += os.system(make_time_dir)

            make_u_dir = f"touch {casePath}/constant/boundaryData/{inletName}/{trunc_time_i}/U;"
            errVal += os.system(make_u_dir)

            make_alpha_dir = f"touch {casePath}/constant/boundaryData/{inletName}/{trunc_time_i}/{alpha_name};"
            errVal += os.system(make_alpha_dir)
        if errVal > 0:
            sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
    except:
        sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
    logger.info("Folder 'boundaryData' was successfully created.")


    # Write the values of flow speed and Volume-Of-Fluid of water per timestep in the newly created boundaryData folder
    # Do not forget to also create boundaryData/inlet/points, containing the coordinates of the inlet points (location where the speed of VOF has to be applied)
    logger.info("Writing inlet definition to folder 'boundaryData'.")
    nPoints = len(coordList[:, 0])
    for i, time_i in enumerate(timeVal):
        trunc_time_i = truncate(time_i, delta_time)

        # Set location of files
        fileLoc_points = f"{casePath}/constant/boundaryData/{inletName}/points"
        U_folder = f"{trunc_time_i}/U"
        fileLoc_U = f"{casePath}/constant/boundaryData/{inletName}/{U_folder}"
        VOFw_folder = f"{trunc_time_i}/{alpha_name}"
        fileLoc_VOFw = f"{casePath}/constant/boundaryData/{inletName}/{VOFw_folder}"

        # Write OpenFOAM-header
        # writeHeader(fileLoc_points, "vectorField", "points")
        # writeHeader(fileLoc_U, "vectorAverageField", "values")
        # writeHeader(fileLoc_VOFw, "scalarAverageField", "values")

        # Write 'points'-file
        f = open(fileLoc_points, 'a+')
        f.write(str(nPoints)+'\n')
        f.write('('+'\n')
        for j in np.arange(nPoints):
            f.write('(' + str(coordList[j, 1]) + ' ' + str(coordList[j, 2]) + ' ' + str(coordList[j, 3]) + ') \n')
        f.write(')'+'\n')
        f.close()
        # Write 'U'-file
        f=open(fileLoc_U,'a+')

        # necessary for vectorAverageField, but not used if 'setAverage'
        # is set to false (checked earlier in this script) TODO
        # f.write('//Average'+'\n')
        # f.write('(0 0 0)'+'\n'+'\n')

        f.write('//Data points'+'\n')
        f.write(str(nPoints)+'\n'+'('+'\n')
        for j in np.arange(nPoints):
            f.write('(' + str(UVal[j, i, 0]) + ' ' + str(UVal[j, i, 1]) + ' ' + str(UVal[j, i, 2]) + ') \n')
        f.write(')'+'\n')
        f.close()
        # Write alpha-file
        f = open(fileLoc_VOFw, 'a+')

        # necessary for vectorAverageField, but not used if 'setAverage'
        # is set to false (checked earlier in this script) TODO
        # f.write('//Average'+'\n')
        # f.write('0'+'\n'+'\n')  #

        f.write('//Data points'+'\n')
        f.write(str(nPoints)+'\n'+'('+'\n')
        for j in np.arange(nPoints):
            f.write(str(VOFwVal[j,i,0])+'\n')
        f.write(')'+'\n')
        f.close()

        # Write OpenFOAM-footers
        writeFooter(fileLoc_points)
        writeFooter(fileLoc_U)
        writeFooter(fileLoc_VOFw)

    # move area postProcess file
    area_path = os.path.join(casePath, "0", "area")
    if os.path.exists(area_path):
        shutil.move(area_path, output_path)


    logger.info("Finished writing boundary condition to folder 'boundaryData'.")
    logger.info("========================End write_bc_openfoam========================")
