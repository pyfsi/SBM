# This script was developed to make a new type of inlet modelling, specifically designed for tube bundle geometries
# but not limited hereto. In this script, (large) bubble shapes are first defined by the user, then randomly selected
# to give a certain mass flux over a unit time and finally written to an appropriate format to serve as a transient
# inlet model for CFD calculations. This script is called automatically by the masterscript 'TubeBundle_master.sh',
# so the user input is channeled to this python script from the bash-script directly.

from utils import np, sys, os, random, PI, NPY_OUT_FOLDER, logging
logger = logging.getLogger(__name__)

class InletModel():
    def __init__(self, UVal, VOFwVal, config):
        self.UVal = UVal
        self.VOFwVal = VOFwVal

        # configuration
        self.startTime = float(config["cfd"]["start_time"])
        self.endTime = float(config["cfd"]["end_time"])
        self.timeStepSize = float(config["cfd"]["delta_time"])
        self.tunit = float(config["sbm"]["t_unit"])
        self.rhog = float(config["cfd"]["rho_g"])
        self.rhol = float(config["cfd"]["rho_l"])
        self.mg_tunit = float(config["sbm"]["mass_gas"])
        self.tol_mg = float(config["sbm"]["mass_gas_tol"])
        self.U = float(config["sbm"]["velocity"])
        # self.mgb_min = float(config["sbm"]["mass_bubble_min"])
        # self.mgb_max = float(config["sbm"]["mass_bubble_max"])
        self.intersectBoundary = config["sbm"]["intersect_boundary"]
        self.intersectBubble = config["sbm"]["intersect_bubble"]

        # store previous bubble centers for next time step
        self.prev_centers = []

    def update(self, UVal, VOFwVal):
        self.UVal = UVal
        self.VOFwVal = VOFwVal

    def define_bubble(self, C_ID, C_t, t, shapeID, mg_bubble, coordList, timeVal, Nshapes, normalInlet):

        startTime = self.startTime
        endTime = self.endTime
        timeStepSize = self.timeStepSize
        tunit = self.tunit
        rhog = self.rhog
        rhol = self.rhol
        U = self.U
        intersectBoundary = self.intersectBoundary
        intersectBubble = self.intersectBubble

        UVal_temp = np.array(self.UVal)
        VOFwVal_temp = np.array(self.VOFwVal)

        C_coord = coordList[C_ID, :]
        C_time = timeVal[C_t]
        if shapeID > (Nshapes-1):
            raise RuntimeError("ShapeID is not within the range 0 to Nshapes-1.")

        # For every bubble, define different requirements for the cell center location and for the scaling factor,
        # which - if not met - causes the end of the current function call. However, always make sure that you have at least
        # one bubble shape that can define a bubble sufficiently small to fall within the set tolerance tol_mg (see further)
        # of the desired amount of gas mg_tunit.
        if shapeID == 0:
            radius_gas = ((3.0*mg_bubble)/(4.0*PI*rhog))**(1.0/3.0)

            # If desired, check that the center point denoted by C_ID and C_t follows a certain set of requirements
            C_checked = True
            timeLoc = C_time - startTime - int((C_time-startTime)/tunit) * tunit  # how many seconds compared to start t_unit

            # Checks below prevents intersection with beginning of t_unit domain
            if timeLoc < radius_gas/U:
                C_checked = False
            if timeLoc > (tunit-radius_gas/U):
                C_checked = False
            if not(C_checked):  # Position of the bubble center (C_ID,C_t) is not OK.
                return False, 0.0

            # Check for each element in the VOFwVal[:, :, 0] whether it's in the bubble to be defined and whether this
            # bubble does not intersect with a previously defined gas bubble. If no old bubble is intersected, change the
            # element VOFwVal and UVal to the appropriate-value; this will be stored in the temporary matrices which will be
            # checked afterwards before updating VOFwVal and UVal
            # Concurrently, integrate the mass of gas you have introduced in the domain.
            mg_checked = 0.0
            mg_bubbleWall = 0.0
            coordCenter = np.array([C_coord[1] - (U * C_time) * normalInlet[0],
                                    C_coord[2] - (U * C_time) * normalInlet[1],
                                    C_coord[3] - (U * C_time) * normalInlet[2]])

            i_mask = np.linalg.norm(coordList[:, 1:4] - C_coord[1:4], axis=1) < radius_gas
            i_list = i_mask.nonzero()[0]
            j_min = int((timeLoc - radius_gas/U)/timeStepSize)
            temp = (timeLoc + radius_gas/U) // timeStepSize
            j_max = int(temp) + 1

            # nested loop could maybe be eliminated by vectorization: boolean arithmetic on the matrices
            n_true_flag = 0
            for i in i_list:
                for j in range(j_min, j_max):
                    time_id = t * int(tunit / timeStepSize) + j
                    coordPoint = np.array(
                        [coordList[i, 1] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[0],
                        coordList[i, 2] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[1],
                        coordList[i, 3] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[2]])
                    if np.linalg.norm(coordPoint-coordCenter) < radius_gas:
                        if self.VOFwVal[i, t * int(
                                tunit / timeStepSize) + j, 0] == 1.0:  # Every cell not yet occupied by bubble
                            VOFwVal_temp[i, t * int(tunit / timeStepSize) + j, 0] = 0.0
                            UVal_temp[i, t * int(tunit / timeStepSize) + j, 0] = U * normalInlet[0]
                            UVal_temp[i, t * int(tunit / timeStepSize) + j, 1] = U * normalInlet[1]
                            UVal_temp[i, t * int(tunit / timeStepSize) + j, 2] = U * normalInlet[2]
                            mg_checked += coordList[i, 4] * U * timeStepSize * rhog
                            mg_bubbleWall = mg_bubbleWall + coordList[i, 4] * U * timeStepSize * rhog
                            n_true_flag += 1
                        elif intersectBubble:
                            mg_bubbleWall = mg_bubbleWall + coordList[
                                i, 4] * U * timeStepSize * rhog  # In this case, a cell was already filled with air, but I will add the mass of air to mg_bubbleWall to be able to check later whether a wall was intersected.
                        else:
                            return False, 0.0

        # Check mass of gas added to the domain: in case intersection with boundary is not allowed (intersectBoundary=False)
        if not(intersectBoundary):
            if mg_bubbleWall < (mg_bubble-np.average(coordList[:, 4])*U*timeStepSize):
                return False, 0.0

        # Save temporary files to permanent files
        self.update(UVal_temp, VOFwVal_temp)

        return True, mg_checked

def inlet_modelling(config):
    logger.info("========================Start inlet_modelling========================")
    print(f"Running inlet_modelling")

    # configuration
    casePath = os.getcwd()
    startTime = float(config["cfd"]["start_time"])
    endTime = float(config["cfd"]["end_time"])
    timeStepSize = float(config["cfd"]["delta_time"])
    tunit = float(config["sbm"]["t_unit"])
    mg_tunit = float(config["sbm"]["mass_gas"])
    tol_mg = float(config["sbm"]["mass_gas_tol"])
    U = float(config["sbm"]["velocity"])
    mgb_min = float(config["sbm"]["mass_bubble_min"])
    mgb_max = float(config["sbm"]["mass_bubble_max"])
    seed = config["sbm"].get("seed", None)
    output_path = os.path.join(casePath, NPY_OUT_FOLDER)

    # set seed for random number generator
    if seed=="None":
        seed = None
    seed_msg = f"Using seed {seed} as {type(seed)} for random number generator"
    print(seed_msg)
    logging.info(seed_msg)
    random.seed(seed)

    if int((endTime-startTime)/tunit) != ((endTime-startTime)/tunit):
        sys.exit("The desired time interval (endTime - startTime) should be a multiple of tunit.")
    if endTime <= startTime:
        sys.exit("The endTime should be larger than the startTime.")
    if (abs(int(tunit/timeStepSize) - tunit/timeStepSize) >= timeStepSize) and (abs((int(tunit/timeStepSize)+1) - tunit/timeStepSize) >= timeStepSize):
        sys.exit("Variable tunit should be a multiple of timeStepSize.")

    # Reading the inlet geometry and normal to the inlet condition
    coordList = np.load(os.path.join(output_path, "inletPython.npy"))
    normalInlet = np.load(os.path.join(output_path, "normalInletPython.npy"))

    # Initializing U and VOFw
    nTimeSteps = int((endTime-startTime)/timeStepSize)+1 # Plus 1 because first time step is included - last time step is not included.
    UVal = np.ones([len(coordList), nTimeSteps, 3])  # initially: "pre-inlet domain" at constant velocity
    for i in np.arange(3):
        UVal[:, :, 0] = U*normalInlet[0]
        UVal[:, :, 1] = U*normalInlet[1]
        UVal[:, :, 2] = U*normalInlet[2]
    VOFwVal = np.ones([len(coordList), nTimeSteps+1, 1])  # "pre-inlet domain" filled with water
    timeVal = np.arange(startTime, endTime+timeStepSize, timeStepSize)  # list of flow times to be defined in this model
    inlet_model = InletModel(UVal, VOFwVal, config)

    # Bubble shapes are defined as 1 function named "define_bubble", comprising a switch based on the shapeID of the bubble
    # (each shape gets its own shapeID and is defined in another part of the switch)
    # Input: centerpoint location cellID in coordList - time instant index in vector timeVal at which cell centre appears -
    # integer 'timeInterval' denoting which time interval is being defined - shapeID-variable: what bubble shape you want to
    # define - amount of mass of gas desired in bubble (scaling factor) - amount of gas in [0,tunit[ which still needs to be
    # defined to get to mg_tunit
    # Output: boolean indicating whether the randomly chosen centerpoint and bubble shape were compatible, if False no
    # bubble will be defined.
    Nshapes = 1  # hard-coded counter used to verify validity of input of "shapeID"-variable - adapt when adding or removing bubble shapes
    probabilityShapes = [1.0] # probability distribution for shape function. sum of all elements must be 0
    if np.sum(probabilityShapes) != 1.0:
        sys.exit('Vector "probabilityShapes" indicating the probability of occurrence of bubble shapes has not been defined correctly.')

    # Selecting bubble shapes based on user input
    # The random generator selects randomly: the bubble shape definition (shapeID) - the center point of a bubble, both in
    # inlet plane and in time (normal direction) - the amount of mass that bubble should have
    # Distribution of the center points is random through entire inlet - restriction to center point location or mass of gas
    # in bubble are defined in the bubble shapes.
    nIntervals = int((endTime-startTime)/tunit)
    logger.info(f"Between startTime {startTime} s and endTime {endTime} s, {nIntervals} intervals of {tunit} s need to be defined.")
    for t in np.arange(nIntervals):
        logger.info(f"Start bubble calculation for time interval {t}")
        iter = 0
        mg_defined = 0.0
        # iterate until mass of inserted gas within tolerance of target mass
        while (mg_tunit-mg_defined) > (tol_mg):
            shapeID = random.randint(0, Nshapes - 1)  # randomly select bubble shape
            # randomly select centerpoint location - determined by cell center ID (2D determined)
            C_ID = random.randint(0, len(coordList) - 1)
            # randomly select centerpoint time location - determined by time step index in timeVal (1D determined)
            C_t = random.randint(t * int(tunit / timeStepSize),
                                 (t + 1) * int(tunit / timeStepSize) - 1)
            # randomly select a scale factor for the bubble
            mg_bubble = random.uniform(min((mgb_min, mg_tunit - mg_defined)),
                                       min((mgb_max, mg_tunit - mg_defined)))

            is_bubble_defined, mg_checked = inlet_model.define_bubble(C_ID, C_t, t,
                                                                      shapeID, mg_bubble, coordList,
                                                                      timeVal, Nshapes, normalInlet)
            if is_bubble_defined:
                mg_defined += mg_checked
                iter = 0
            else:
                iter = iter+1
                residual = mg_tunit-mg_defined
                # print(f"\t mg_bubble {mg_bubble:.10f} \t Residual {residual:.10f}")
            if iter > 1000:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations. Bubble generation with the current configuration file not possible.")
        logger.info(f"\t Mass of inserted gas: {mg_defined} kg. (Target mass: {mg_tunit} kg).")


    logger.info(f"Inlet model iteration loop ended.")

    # Writing the profile to be used in OpenFOAM
    logger.info("Saving inlet profile to Python (numpy) npy-files.")
    # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of velocity
    np.save(os.path.join(output_path, "inletDefinition-U.npy"), inlet_model.UVal)
    # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of VOFw
    np.save(os.path.join(output_path, "inletDefinition-VOFw.npy"), inlet_model.VOFwVal)
    # List containing the time instants where U and alpha.water are defined
    np.save(os.path.join(output_path, "inletDefinition-time.npy"), timeVal)
    logger.info("Inlet profile saved in Python (numpy) npy-files.")

    # Check: convert to file compatible with ParaView to visualize your pre-inlet domain.
    logger.info("Saving inlet profile to CSV-files. ")
    files = [os.path.join(output_path, "inletDefinition-VOFw.csv"),
             os.path.join(output_path, "inletDefinition-Ux.csv"),
             os.path.join(output_path, "inletDefinition-Uy.csv"),
             os.path.join(output_path, "inletDefinition-Uz.csv"),
            ]
    toWrite = [inlet_model.VOFwVal[:, :, 0], inlet_model.UVal[:, :, 0],
               inlet_model.UVal[:, :, 1], inlet_model.UVal[:, :, 2]]
    for fi in np.arange(len(files)):
        f = open(files[fi], 'w')
        f.write('x coord,y coord,z coord,value \n')
        for i in np.arange(len(coordList[:, 0])):
            for j in np.arange(len(timeVal)):
                coordPoint = np.array([coordList[i, 1] - (U * timeVal[j]) * normalInlet[0],
                                    coordList[i, 2] - (U * timeVal[j]) * normalInlet[1],
                                    coordList[i, 3] - (U * timeVal[j]) * normalInlet[2]])
                f.write(str(coordPoint[0]) + ',' + str(coordPoint[1]) + ',' + str(coordPoint[2]) + ',' + str(
                    toWrite[fi][i, j]) + '\n')
        f.close()
    logger.info("Inlet profile saved to CSV-files.")

    logger.info("========================End inlet_modelling========================")
