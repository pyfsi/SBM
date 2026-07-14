# This script is developed to read in an OpenFOAM-case, search for a specific inlet boundary and rewrite the cell
# centres and faces on that boundary to an array to be used for the inlet modelling in a later stage (not in this
# script). This script is called automatically by the masterscript 'TubeBundle_master.sh', so the user input is
# channeled to this python script from the bash-script directly.

from utils import np, os, sys, linecache, NPY_OUT_FOLDER, logging
logger = logging.getLogger(__name__)

def read_inlet_openfoam(config):
    logger.info("========================Start read_inlet_openfoam========================")
    print(f"Running read_inlet_openfoam")

    # alias for config variables
    casePath = os.getcwd()
    moduleCFD = str(config["packages"]["cfd_program"] + "/" + config["packages"]["cfd_version"])
    openfoam_type = str(config.get("openfoam_type"))
    dimesion = int(config["cfd"]["dimension"])
    startTime = str(config["cfd"]["start_time"])
    inletName = str(config["cfd"]["inlet_name"])

    # create directory for sbm output if it doesn't exist
    output_path = os.path.join(casePath, NPY_OUT_FOLDER)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Read inlet boundary from OpenFOAM
    logger.info("Starting 'postProcess' module of OpenFOAM")
    if openfoam_type=="com":
        os.system("cd " + casePath + "; ml " + moduleCFD + "; source $FOAM_BASH; postProcess -time " + startTime + ";")
    if openfoam_type=="org":
        os.system("cd " + casePath + "; ml " + moduleCFD + "; source $FOAM_BASH; foamPostProcess; foamPostProcess -func writeCellCentres;")
    logger.info("Finished 'postProcess'")

    # filenames
    if openfoam_type=="org":
        source_filenames = ["Ccx", "Ccy", "Ccz", "area"]
    else:
        source_filenames = ["Cx", "Cy", "Cz", "area"]
    logger.info("Loading inlet cell coordinates and face areas into Python.")
    for i, filename in enumerate(source_filenames):
        sourceFile = f"{casePath}/{startTime}/{filename}"

        try:
            # find line where inlet-list starts; you should redo this for all coordinates as 'writeCellCentres' - when
            # possible - replaces lists of uniform values by a single uniform value statement. If the latter is the case,
            # the number of lines you should have in a complete list is obtained from the number of list elements read in
            # another coordinate-file (it's not possible to have >1 'uniform value' coordinate files). The incomplete
            # uniform value-files are subsequently (outside the try-loop) converted to a list of the same size as the
            # non-uniform list coordinates.
            if not os.path.exists(sourceFile):
                sys.exit(f"{sourceFile} not found")

            os.system("cd " + casePath + "; grep -nr " + inletName + " " + sourceFile + " | cut -d : -f 1 > lineNr")
            lineNameNr = int(open(casePath + "/lineNr", 'r').readline())
            lineStartNr = lineNameNr + 6  # In case of non-uniform list, this is where the list of values in the sourceFile starts
            rowsNrIndex = lineNameNr + 4  # On this line, the number of cell centers on the inlet is stated
            os.system("cd " + casePath + "; awk NR==" + str(rowsNrIndex) + " " + sourceFile + " > rowsNr")
            rowsNr = int(open(casePath + "/rowsNr", 'r').readline())
            os.system("cd " + casePath + "; rm lineNr rowsNr")
            tempCoordFile = np.ones([rowsNr, 1]) * float("inf")
            for j in np.arange(rowsNr):
                tempCoordFile[j, 0] = float(linecache.getline(sourceFile, lineStartNr + j))
        except ValueError: # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking
            if 'rowsNr' not in locals(): # If the first coordinate-file has a uniform value, the variable 'rowsNr' does not exist, so you should check whether this variable exists
                checkFile = casePath + "/" + startTime + "/ccy" # if 'rowsNr' does not exist, read the second sourceFile to know the number of rows
                os.system("cd " + casePath + "; grep -nr " + inletName + " " + checkFile + " | cut -d : -f 1 > lineNr")
                lineNameNr_CF = int(open(casePath+"/lineNr",'r').readline())
                rowsNrIndex_CF = lineNameNr_CF+4 # On this line, the number of cell centers on the inlet is stated
                os.system("cd " + casePath + "; awk NR==" + str(rowsNrIndex_CF) + " " + checkFile + " > rowsNr")
                rowsNr = int(open(casePath+"/rowsNr", 'r').readline())
                os.system("cd " + casePath + "; rm lineNr rowsNr")
            indexUV = lineNameNr + 3
            os.system("cd " + casePath + "; awk NR==" + str(indexUV) + " " + sourceFile + " > unifValue")
            unifValue = float(open(casePath+"/unifValue", 'r').readline().split()[-1][0:-1]) # First '-1'm akes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
            os.system("cd " + casePath + "; rm unifValue")
            tempCoordFile = np.ones([rowsNr, 1])*float("inf")
            for j in np.arange(rowsNr):
                tempCoordFile[j, 0] = unifValue

        if i == 0:
            coordList = np.ones([rowsNr, 5]) * float("inf")  # ID - X - Y - Z - Area
            coordList[:, 0] = np.arange(rowsNr)

        coordList[:, (i + 1)] = tempCoordFile[:, 0]


    # Check that all values are inserted correctly
    for i in np.arange(rowsNr):
        for j in np.arange(5):
            if coordList[i, j] > 1e15:
                sys.exit("Not all values are correctly read into the Python-script 'TubeBundle_ReadInlet'")
    logger.info("Completed loading of cell coordinates and face areas into Python.")

    # Determine the normal pointing of the inlet, pointing INTO the domain
    # For this purpose, find three linearly independent points in coordList
    logger.info("Calculating inlet normals.")
    if dimesion == 3:
        tol_product = 0.01
        point1 = coordList[0, 1:4]
        point2 = coordList[1, 1:4]
        i = 2
        point3 = coordList[i, 1:4]
        enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
        denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))
        while (enum / denom) < tol_product:
            i = i + 1
            point3 = coordList[i, 1:4]

            enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
            denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))

        normal_vec = np.cross(point2 - point1, point3 - point1)
        normalInlet = normal_vec / np.linalg.norm(normal_vec)

        # Need one more point from the domain to determine the correct orientation of the inlet normal
        os.system("cd " + casePath + "; grep -nr '(' constant/polyMesh/points | head -n 1 | cut -d : -f 1 > lineNr")
        lineNameNr = int(open(casePath+"/lineNr", 'r').readline())
        lineNr = lineNameNr+1  # First point that is defined
        os.system("cd " + casePath + "; rm lineNr ")
        i = 0
        f = open(casePath+"/constant/polyMesh/points", 'r')
        while i < lineNr:
            i = i+1
            f.readline()
        pointDomain = np.double(np.array(f.readline()[1:-2].split(" ")))
        while np.linalg.norm(np.dot(pointDomain-point1, normalInlet))/(np.linalg.norm(pointDomain-point1)) < tol_product: #find pointDomain which is not in the inlet plane - normalInlet already has norm equal to 1
            pointDomain = np.double(np.array(f.readline()[1:-2].split(" ")))
        if np.dot(pointDomain-point1, normalInlet) < 0:
            normalInlet = (-1.0)*normalInlet
        f.close()
    elif dimesion == 2:
        normalInlet = np.zeros([3])
        print("The normal to the inlet cannot be calculated directly. Please give the x-, y- and z-coordinates of the normal vector in the following prompts.")
        axis = ["x", "y", "z"]
        for i, ax_i in enumerate(axis):
            coordOK = False
            while not coordOK:
                try:
                    print(f"Please provide the {ax_i}-component of the normal vector: ")
                    temp_coord = float(input())
                    coordOK = True
                except ValueError:
                    print("Please provide a float value. Try again.")
            normalInlet[i] = temp_coord
        normalInlet = (1/np.linalg.norm(normalInlet))*normalInlet
    else:
        sys.exit("Number of dimesion should be either 2 or 3.")

    logger.info(f"Completed calculating normal to the inlet, pointing into the domain: {normalInlet}")

    # Save inlet and normal in Python Numpy-array format
    np.save(os.path.join(output_path, "inletPython.npy"), coordList)
    np.save(os.path.join(output_path, "normalInletPython.npy"), normalInlet)

    logger.info("========================End read_inlet_openfoam========================")