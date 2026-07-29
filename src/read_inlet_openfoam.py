from utils import np, os, sys, linecache, SBM_OUTPUT, logging
logger = logging.getLogger(__name__)

TOLERANCE = 0.01

def read_inlet_openfoam(config):
    logger.info("========================Start read_inlet_openfoam========================")
    print(f"Running read_inlet_openfoam")

    # alias for config variables
    case_path = os.getcwd()
    cfd_module = str(config["packages"]["cfd_program"] + "/" + config["packages"]["cfd_version"])
    openfoam_type = str(config.get("openfoam_type"))
    dimesion = int(config["cfd"]["dimension"])
    start_time = str(config["cfd"]["start_time"])
    inlet_name = str(config["cfd"]["inlet_name"])

    # create directory for sbm output if it doesn't exist
    output_path = os.path.join(case_path, SBM_OUTPUT)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Read inlet boundary from OpenFOAM
    logger.info("Starting 'postProcess' module of OpenFOAM")
    if openfoam_type=="com":
        os.system("cd " + case_path + "; ml " + cfd_module + "; source $FOAM_BASH; postProcess -time " + start_time + ";")
    if openfoam_type=="org":
        os.system("cd " + case_path + "; ml " + cfd_module + "; source $FOAM_BASH; foamPostProcess; foamPostProcess -func writeCellCentres;")
    logger.info("Finished 'postProcess'")

    # filenames
    if openfoam_type=="org":
        inlet_geometry_filenames = ["Ccx", "Ccy", "Ccz", "area"]
    else:
        inlet_geometry_filenames = ["Cx", "Cy", "Cz", "area"]
    logger.info("Loading inlet cell coordinates and face areas into Python.")
    for i, file_i in enumerate(inlet_geometry_filenames):
        source_file = f"{case_path}/{start_time}/{file_i}"

        try:
            # find line where inlet-list starts; you should redo this for all coordinates as 'writeCellCentres' - when
            # possible - replaces lists of uniform values by a single uniform value statement. If the latter is the case,
            # the number of lines you should have in a complete list is obtained from the number of list elements read in
            # another coordinate-file (it's not possible to have >1 'uniform value' coordinate files). The incomplete
            # uniform value-files are subsequently (outside the try-loop) converted to a list of the same size as the
            # non-uniform list coordinates.
            if not os.path.exists(source_file):
                sys.exit(f"{source_file} not found")

            os.system("cd " + case_path + "; grep -nr " + inlet_name + " " + source_file + " | cut -d : -f 1 > lineNr")
            lineNameNr = int(open(case_path + "/lineNr", 'r').readline())
            lineStartNr = lineNameNr + 6  # In case of non-uniform list, this is where the list of values in the source_file starts
            rowsNrIndex = lineNameNr + 4  # On this line, the number of cell centers on the inlet is stated
            os.system("cd " + case_path + "; awk NR==" + str(rowsNrIndex) + " " + source_file + " > rowsNr")
            rowsNr = int(open(case_path + "/rowsNr", 'r').readline())
            os.system("cd " + case_path + "; rm lineNr rowsNr")
            tempCoordFile = np.ones([rowsNr, 1]) * float("inf")
            for j in np.arange(rowsNr):
                tempCoordFile[j, 0] = float(linecache.getline(source_file, lineStartNr + j))
        except ValueError: # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking
            if 'rowsNr' not in locals(): # If the first coordinate-file has a uniform value, the variable 'rowsNr' does not exist, so you should check whether this variable exists
                checkFile = case_path + "/" + start_time + "/ccy" # if 'rowsNr' does not exist, read the second source_file to know the number of rows
                os.system("cd " + case_path + "; grep -nr " + inlet_name + " " + checkFile + " | cut -d : -f 1 > lineNr")
                lineNameNr_CF = int(open(case_path+"/lineNr",'r').readline())
                rowsNrIndex_CF = lineNameNr_CF+4 # On this line, the number of cell centers on the inlet is stated
                os.system("cd " + case_path + "; awk NR==" + str(rowsNrIndex_CF) + " " + checkFile + " > rowsNr")
                rowsNr = int(open(case_path+"/rowsNr", 'r').readline())
                os.system("cd " + case_path + "; rm lineNr rowsNr")
            indexUV = lineNameNr + 3
            os.system("cd " + case_path + "; awk NR==" + str(indexUV) + " " + source_file + " > unifValue")
            unifValue = float(open(case_path+"/unifValue", 'r').readline().split()[-1][0:-1]) # First '-1'm akes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
            os.system("cd " + case_path + "; rm unifValue")
            tempCoordFile = np.ones([rowsNr, 1])*float("inf")
            for j in np.arange(rowsNr):
                tempCoordFile[j, 0] = unifValue

        if i == 0:
            coord_list = np.ones([rowsNr, 5]) * float("inf")  # ID - X - Y - Z - Area
            coord_list[:, 0] = np.arange(rowsNr)

        coord_list[:, (i + 1)] = tempCoordFile[:, 0]


    # Check that all values are inserted correctly
    for i in np.arange(rowsNr):
        for j in np.arange(5):
            if coord_list[i, j] > 1e15:
                sys.exit("Not all values are correctly read into the Python-script 'TubeBundle_ReadInlet'")
    logger.info("Completed loading of cell coordinates and face areas into Python.")

    # Calculate inlet normal
    logger.info("Calculating inlet normals.")
    if dimesion == 3:
        point1 = coord_list[0, 1:4]
        point2 = coord_list[1, 1:4]
        i = 2
        point3 = coord_list[i, 1:4]
        enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
        denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))
        while (enum / denom) < TOLERANCE:
            i = i + 1
            point3 = coord_list[i, 1:4]
            enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
            denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))

        normal_vec = np.cross(point2 - point1, point3 - point1)
        normal_inlet = normal_vec / np.linalg.norm(normal_vec)

        # Need one more point from the domain to determine the correct orientation of the inlet normal
        os.system("cd " + case_path + "; grep -nr '(' constant/polyMesh/points | head -n 1 | cut -d : -f 1 > lineNr")
        lineNameNr = int(open(case_path+"/lineNr", 'r').readline())
        lineNr = lineNameNr+1  # First point that is defined
        os.system("cd " + case_path + "; rm lineNr ")
        i = 0
        points_path = os.path.join(case_path, "constant", "polyMesh", "points")
        with open(points_path, 'r') as f:
            while i < lineNr:
                i = i+1
                f.readline()
            point_inside_domain = np.double(np.array(f.readline()[1:-2].split(" ")))

            #find point_inside_domain which is not in the inlet plane
            while abs(np.dot(point_inside_domain-point1, normal_inlet))/(np.linalg.norm(point_inside_domain-point1)) < TOLERANCE:
                point_inside_domain = np.double(np.array(f.readline()[1:-2].split(" ")))

            # switch sign if normal_inlet pointing outside
            if np.dot(point_inside_domain-point1, normal_inlet) < 0:
                normal_inlet = (-1.0)*normal_inlet
    elif dimesion == 2:
        normal_inlet = np.zeros([3])
        print("The normal to the inlet cannot be calculated directly. Please input the x-, y- and z-coordinates of the normal vector.")
        axis = ["x", "y", "z"]
        for i, ax_i in enumerate(axis):
            is_normal_correct = False
            while not is_normal_correct:
                try:
                    print(f"Please provide the {ax_i}-component of the normal vector: ")
                    temp_normal = float(input())
                    is_normal_correct = True
                except ValueError:
                    print("Input values is not expected type. Please provide a float.")
            normal_inlet[i] = temp_normal
        normal_inlet = (1/np.linalg.norm(normal_inlet))*normal_inlet
    else:
        sys.exit("Number of dimesion should be either 2 or 3.")

    logger.info(f"Completed calculating normal to the inlet, pointing into the domain: {normal_inlet}")

    # Save inlet and normal in Python Numpy-array format
    np.save(os.path.join(output_path, "inletPython.npy"), coord_list)
    np.save(os.path.join(output_path, "normalInletPython.npy"), normal_inlet)

    logger.info("========================End read_inlet_openfoam========================")