# This script was developed to make a new type of inlet modelling, specifically designed for tube bundle geometries
# but not limited hereto. In this script, (large) bubble shapes are first defined by the user, then randomly selected
# to give a certain mass flux over a unit time and finally written to an appropriate format to serve as a transient
# inlet model for CFD calculations. This script is called automatically by the masterscript 'TubeBundle_master.sh',
# so the user input is channeled to this python script from the bash-script directly.


# import of utilities
import math
import numpy as np
import sys
import random  # package with random generator
import yaml

print("Starting inlet modelling script.\n")

# Read configuration file
with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.SafeLoader)
casePath = config["case_path"]
startTime = float(config["simulation_parameters"]["start_time"])
endTime = float(config["simulation_parameters"]["end_time"])
timeStepSize = float(config["simulation_parameters"]["delta_time"])
tunit = float(config["sbm_parameters"]["t_unit"])
inletName = config["simulation_parameters"]["inlet_name"]
rhog = float(config["simulation_parameters"]["rho_g"])
rhol = float(config["simulation_parameters"]["rho_l"])
mg_tunit = float(config["sbm_parameters"]["mass_bubble"])
tol_mg = float(config["sbm_parameters"]["mass_bubble_tol"])
U = float(config["sbm_parameters"]["velocity"])
mgb_min = float(config["sbm_parameters"]["mass_bubble_min"])
mgb_max = float(config["sbm_parameters"]["mass_bubble_max"])
intersectBoundary = config["sbm_parameters"]["intersect_boundary"]
intersectBubble = config["sbm_parameters"]["intersect_bubble"]

if int((endTime-startTime)/tunit) != ((endTime-startTime)/tunit):
    sys.exit("The desired time interval (endTime - startTime) should be a multiple of tunit.")
if endTime <= startTime:
    sys.exit("The endTime should be larger than the startTime.")
if (abs(int(tunit/timeStepSize) - tunit/timeStepSize) >= timeStepSize) and (abs((int(tunit/timeStepSize)+1) - tunit/timeStepSize) >= timeStepSize):
    sys.exit("Variable tunit should be a multiple of timeStepSize.")  

# Reading the inlet geometry and normal to the inlet condition prepared with 'TubeBundle_readInlet_<CFD-programme>.py' 
coordList = np.load(casePath + "/inletPython.npy")
normalInlet = np.load(casePath + "/normalInletPython.npy")

# Initializing U and VOFw
nTimeSteps = int((endTime-startTime)/timeStepSize)+1 # Plus 1 because first time step is included - last time step is not included.
UVal = np.ones([len(coordList), nTimeSteps, 3])  # initially: "pre-inlet domain" at constant velocity
for i in np.arange(3):
    UVal[:, :, 0] = U*normalInlet[0]
    UVal[:, :, 1] = U*normalInlet[1]
    UVal[:, :, 2] = U*normalInlet[2]
VOFwVal = np.ones([len(coordList), nTimeSteps, 1])  # initially: "pre-inlet domain" filled with water
timeVal = np.arange(startTime, endTime, timeStepSize)  # list of flow times to be defined in this model


# Creating bubble shapes - under the hood, so hard-coded shapes
# Bubble shapes are defined as 1 function named "bubbleShape", comprising a switch based on the shapeID of the bubble
# (each shape gets its own shapeID and is defined in another part of the switch)
# Input: centerpoint location cellID in coordList - time instant index in vector timeVal at which cell centre appears -
# integer 'timeInterval' denoting which time interval is being defined - shapeID-variable: what bubble shape you want to
# define - amount of mass of gas desired in bubble (scaling factor) - amount of gas in [0,tunit[ which still needs to be
# defined to get to mg_tunit
# Output: boolean indicating whether the randomly chosen centerpoint and bubble shape were compatible, if False no
# bubble will be defined.
Nshapes = 1  # hard-coded counter used to verify validity of input of "shapeID"-variable - adapt when adding or removing bubble shapes
probabilityShapes = [1.0] # hard-coded probability distribution of the bubble shapes - el. 0 = probability of bubble shape 0 .... // It is checked that the sum of elements equals zero.
if np.sum(probabilityShapes) != 1.0:
    sys.exit('Vector "probabilityShapes" indicating the probability of occurrence of bubble shapes has not been defined correctly.')

def bubbleShape(C_ID, C_t, timeInterval, shapeID, mgb):
    # define UVal and VOFwVal to be global such that these matrices can be altered directly by this function 
    # As coordList will not be adapted in this function, 
    # it does not need to be defined as global (Python automatically looks for coordList definition outside of function)
    global UVal, VOFwVal
     
    UVal_temp = np.array(UVal)
    VOFwVal_temp = np.array(VOFwVal)
    C_coord = coordList[C_ID, :]
    C_time = timeVal[C_t]
    if shapeID > (Nshapes-1):
        sys.exit("Fatal error. Shape generator asks for non-existing bubble shape (Nshapes is too low or too few shapes have been defined).")
       
    # Start of switch: every switch input is another bubble shape
    # For every bubble, you can define different requirements for the cell center location and for the scaling factor,
    # which - if not met - causes the end of the current function call. However, always make sure that you have at least
    # one bubble shape that can define a bubble sufficiently small to fall within the set tolerance tol_mg (see further)
    # of the desired amount of gas mg_tunit.
    if shapeID == 0:
        rg = ((3.0*mgb)/(4.0*math.pi*rhog))**(1.0/3.0)

        # If desired, check that the center point denoted by C_ID and C_t follows a certain set of requirements
        C_checked = True
        timeLoc = C_time - startTime - int((C_time-startTime)/tunit) * tunit  # how many seconds compared to start t_unit

        # Checks below prevents intersection with beginning of t_unit domain
        if timeLoc < rg/U:
            C_checked = False
        if timeLoc > (tunit-rg/U):
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

        i_mask = np.linalg.norm(coordList[:, 1:4] - C_coord[1:4], axis=1) < rg
        i_list = i_mask.nonzero()[0]

        j_min = int((timeLoc - rg/U)/timeStepSize)
        j_max = int(math.ceil((timeLoc + rg/U)/timeStepSize)) + 1
        
        # nested loop could maybe be eliminated by vectorization: boolean arithmetic on the matrices
        n_true_flag = 0
        for i in i_list:
            for j in range(j_min, j_max):
                coordPoint = np.array(
                    [coordList[i, 1] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[0],
                     coordList[i, 2] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[1],
                     coordList[i, 3] - (U * timeVal[t * int(tunit / timeStepSize) + j]) * normalInlet[2]])
                if np.linalg.norm(coordPoint-coordCenter) < rg:
                    if VOFwVal[i, timeInterval * int(
                            tunit / timeStepSize) + j, 0] == 1.0:  # Every cell not yet occupied by bubble
                        VOFwVal_temp[i, timeInterval * int(tunit / timeStepSize) + j, 0] = 0.0
                        UVal_temp[i, timeInterval * int(tunit / timeStepSize) + j, 0] = U * normalInlet[0]
                        UVal_temp[i, timeInterval * int(tunit / timeStepSize) + j, 1] = U * normalInlet[1]
                        UVal_temp[i, timeInterval * int(tunit / timeStepSize) + j, 2] = U * normalInlet[2]
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
        if mg_bubbleWall < (mgb-np.average(coordList[:, 4])*U*timeStepSize):
            return False, 0.0
#     if TrialBoolean:
#         print("mg_checked: "+str(mg_checked))
#         print("mgb: "+str(mgb))

    # Save temporary files to permanent files
    UVal = UVal_temp
    VOFwVal = VOFwVal_temp

    return True, mg_checked
           
           
# Selecting bubble shapes based on user input
# The random generator selects randomly: the bubble shape definition (shapeID) - the center point of a bubble, both in
# inlet plane and in time (normal direction) - the amount of mass that bubble should have
# Distribution of the center points is random through entire inlet - restriction to center point location or mass of gas
# in bubble are defined in the bubble shapes.
# The former is 
nIntervals = int((endTime-startTime)/tunit)  # Number of intervals [0,tunit[
print("Between startTime " + str(startTime) + "s and endTime " + str(endTime) + "s, " + str(
    nIntervals) + " intervals of " + str(tunit) + "s need to be defined.")
for t in np.arange(nIntervals):
    print(f"Start bubble calculation for time interval {t}")
    iter = 0
    mg_defined = 0.0  # Variable checking the amount of gas already defined
    while (mg_tunit-mg_defined) > (tol_mg):
        shapeID = random.randint(0, Nshapes - 1)  # randomly select bubble shape
        # randomly select centerpoint location - determined by cell center ID (2D determined)
        C_ID = random.randint(0, len(coordList) - 1)
        # randomly select centerpoint time location - determined by time step index in timeVal (1D determined)
        C_t = random.randint(t * int(tunit / timeStepSize), (t + 1) * int(tunit / timeStepSize) - 1)
        # randomly select a scale factor for the bubble you are creating
        mg_bubble = random.uniform(min((mgb_min, mg_tunit - mg_defined)), min((mgb_max, mg_tunit - mg_defined)))

        bubbleDefined, mg_checked = bubbleShape(C_ID, C_t, t, shapeID, mg_bubble)
        if bubbleDefined:
            mg_defined += mg_checked
            iter = 0
        else:
            iter = iter+1
            temp = mg_tunit-mg_defined
            print(f"\t mg_bubble {mg_bubble:.10f} \t Residual {temp:.10f}")
        if iter > 1000:
            sys.exit("Forced exit: 1000 trials of bubble definition have failed; system is ill-defined.")
    print("\t => Time interval " + str(t) + " has been defined: "+str(mg_defined)+"kg of gas was inserted. (desired: "+str(mg_tunit)+"kg).")
print("Inlet was modelled successfully. \n")     

# Writing the profile to be used in OpenFOAM
print("Saving inlet profile to Python (numpy) npy-files.")
np.save(casePath+'/inletDefinition-U.npy', UVal) # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of velocity
np.save(casePath+'/inletDefinition-VOFw.npy', VOFwVal) # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of VOFw
np.save(casePath+'/inletDefinition-time.npy', timeVal) # List containing the time instants where U and alpha.water are defined
print("Inlet profile saved in Python (numpy) npy-files. \n")

# Check: convert to file compatible with ParaView to visualize your pre-inlet domain.
print("Saving inlet profile to CSV-files. ")
files = [casePath+'/inletDefinition-VOFw.csv', casePath+'/inletDefinition-Ux.csv', casePath+'/inletDefinition-Uy.csv', casePath+'/inletDefinition-Uz.csv']
toWrite = [VOFwVal[:, :, 0], UVal[:, :, 0], UVal[:, :, 1], UVal[:, :, 2]]
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
print("Inlet profile saved to CSV-files. ")


print("Script 'inletModelling' completed. \n")
