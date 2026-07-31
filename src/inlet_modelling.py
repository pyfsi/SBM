from utils import np, sys, os, random, PI, SBM_OUTPUT, logging
logger = logging.getLogger(__name__)

# constants
N_SHAPES = 1

class InletModel():
    def __init__(self, config, coord_list, normal_inlet):

        # configuration
        self.t_start = float(config["cfd"]["start_time"])
        self.t_end = float(config["cfd"]["end_time"])
        self.dt = float(config["cfd"]["delta_time"])
        self.dt_insert = float(config["sbm"]["t_unit"])
        self.rhog = float(config["cfd"]["rho_g"])
        self.rhol = float(config["cfd"]["rho_l"])
        self.mg_per_insert = float(config["sbm"]["mass_gas"])
        self.tol_mg = float(config["sbm"]["mass_gas_tol"])
        self.velocity_bc = float(config["sbm"]["velocity"])
        self.mgb_min = float(config["sbm"]["mass_bubble_min"])
        self.mgb_max = float(config["sbm"]["mass_bubble_max"])
        self.intersect_boundary = config["sbm"]["intersect_boundary"]
        self.intersect_bubble = config["sbm"]["intersect_bubble"]

        # inlet geometry information
        self.coord_list = coord_list
        self.normal_inlet = normal_inlet

        # store previous bubble centers for next time step
        self.prev_centers = []

        # Initializing velocity and phase fraction arrays
        n_time_step = int((self.t_end - self.t_start)/self.dt)+1
        self.velocity = np.ones([len(coord_list), n_time_step, 3])
        self.velocity[:, :, 0] = self.velocity_bc * normal_inlet[0]
        self.velocity[:, :, 1] = self.velocity_bc * normal_inlet[1]
        self.velocity[:, :, 2] = self.velocity_bc * normal_inlet[2]
        self.alpha = np.ones([len(coord_list), n_time_step+1, 1])
        self.time = np.arange(self.t_start, self.t_end+self.dt, self.dt)

    def update(self, velocity, alpha):
        self.velocity = velocity
        self.alpha = alpha

    def define_bubble(self, C_ID, C_t, t, shapeID, mg_bubble):

        t_start = self.t_start
        t_end = self.t_end
        dt = self.dt
        dt_insert = self.dt_insert
        rhog = self.rhog
        rhol = self.rhol
        velocity_bc = self.velocity_bc
        intersect_boundary = self.intersect_boundary
        intersect_bubble = self.intersect_bubble
        coord_list = self.coord_list
        normal_inlet = self.normal_inlet

        # temporary storage arrays
        velocity_temp = self.velocity
        alpha_temp = self.alpha

        C_coord = coord_list[C_ID, :]
        C_time = self.time[C_t]
        if shapeID > (N_SHAPES-1):
            raise RuntimeError("ShapeID is not within the range 0 to N_SHAPES-1.")

        if shapeID == 0:
            radius_gas = ((3.0*mg_bubble)/(4.0*PI*rhog))**(1.0/3.0)

            # Check that the center point denoted by C_ID and C_t follows a certain set of requirements
            C_checked = True
            time_loc = C_time - t_start - int((C_time-t_start)/dt_insert) * dt_insert

            # Checks below prevents intersection with beginning of t_unit domain
            if time_loc < radius_gas/velocity_bc:
                C_checked = False
            if time_loc > (dt_insert-radius_gas/velocity_bc):
                C_checked = False
            if not C_checked:
                return False, 0.0

            mg_defined = 0.0
            mg_bubble_wall = 0.0
            bubble_center = np.array([C_coord[1] - (velocity_bc * C_time) * normal_inlet[0],
                                    C_coord[2] - (velocity_bc * C_time) * normal_inlet[1],
                                    C_coord[3] - (velocity_bc * C_time) * normal_inlet[2]])

            # get space and time index (i and j) of bounding cylinder
            i_mask = np.linalg.norm(coord_list[:, 1:4] - C_coord[1:4], axis=1) < radius_gas
            i_list = i_mask.nonzero()[0]
            j_min = int((time_loc - radius_gas/velocity_bc)/dt)
            temp = (time_loc + radius_gas/velocity_bc) // dt
            j_max = int(temp) + 1

            n_true_flag = 0
            for i in i_list:
                for j in range(j_min, j_max):
                    time_id = t * int(dt_insert / dt) + j
                    coordPoint = np.array(
                        [coord_list[i, 1] - (velocity_bc * self.time[time_id]) * normal_inlet[0],
                        coord_list[i, 2] - (velocity_bc * self.time[time_id]) * normal_inlet[1],
                        coord_list[i, 3] - (velocity_bc * self.time[time_id]) * normal_inlet[2]])
                    if np.linalg.norm(coordPoint-bubble_center) < radius_gas:
                        if self.alpha[i, time_id, 0] == 1.0:
                            alpha_temp[i, time_id, 0] = 0.0
                            velocity_temp[i, time_id, :] = velocity_bc * normal_inlet[:]

                            mg_defined += coord_list[i, 4] * velocity_bc * dt * rhog
                            mg_bubble_wall += coord_list[i, 4] * velocity_bc * dt * rhog
                            n_true_flag += 1
                        elif intersect_bubble:
                            mg_bubble_wall = mg_bubble_wall + coord_list[
                                i, 4] * velocity_bc * dt * rhog
                        else:
                            return False, 0.0

        if not(intersect_boundary):
            if mg_bubble_wall < (mg_bubble-np.average(coord_list[:, 4])*velocity_bc*dt):
                return False, 0.0

        # Save temporary files to permanent files
        self.update(velocity_temp, alpha_temp)

        return True, mg_defined

    def insert_bubbles(self, t_insert_idx) -> float:
        # aliases
        mg_per_insert = self.mg_per_insert
        tol_mg = self.tol_mg
        coord_list = self.coord_list
        dt_insert = self.dt_insert
        dt = self.dt
        mgb_min = self.mgb_min
        mgb_max = self.mgb_max

        # clamp time insert
        n_dt_per_insert = int(dt_insert / dt)
        C_t_min = t_insert_idx * n_dt_per_insert
        C_t_max = (t_insert_idx + 1) * n_dt_per_insert - 1

        iter = 0
        mg_inserted = 0.0
        # iterate until mass of inserted gas within tolerance of target mass
        while abs(mg_per_insert-mg_inserted) > (tol_mg):
            shapeID = 0 #random.randint(0, N_SHAPES - 1)  # randomly select bubble shape
            # sample centerpoint location
            C_ID = random.randint(0, len(coord_list) - 1)
            # sample centerpoint time location
            C_t = random.randint(t_insert_idx * int(dt_insert / dt),
                                    (t_insert_idx + 1) * int(dt_insert / dt) - 1)
            # sample bubble mass
            mg_bubble = random.uniform(min((mgb_min, mg_per_insert - mg_inserted)),
                                        min((mgb_max, mg_per_insert - mg_inserted)))

            is_bubble_defined, mg_defined = self.define_bubble(C_ID, C_t, t_insert_idx, shapeID, mg_bubble)
            if is_bubble_defined:
                mg_inserted += mg_defined
                iter = 0
            else:
                iter = iter+1
                residual = mg_per_insert-mg_inserted
                # print(f"\t mg_bubble {mg_bubble:.10f} \t Residual {residual:.10f}")
            if iter > 1000:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations.")

        return mg_inserted

def inlet_modelling(config):
    logger.info("========================Start inlet_modelling========================")
    print(f"Running inlet_modelling")

    # configuration
    casePath = os.getcwd()
    t_start = float(config["cfd"]["start_time"])
    t_end = float(config["cfd"]["end_time"])
    dt = float(config["cfd"]["delta_time"])
    dt_insert = float(config["sbm"]["t_unit"])
    mg_per_insert = float(config["sbm"]["mass_gas"])
    velocity_bc = float(config["sbm"]["velocity"])
    seed = config["sbm"].get("seed", None)
    output_path = os.path.join(casePath, SBM_OUTPUT)

    # set seed for random number generator
    if seed=="None":
        seed = None
    seed_msg = f"Using seed {seed} as {type(seed)} for random number generator"
    print(seed_msg)
    logging.info(seed_msg)
    random.seed(seed)

    if int((t_end-t_start)/dt_insert) != ((t_end-t_start)/dt_insert):
        sys.exit("The desired time interval (t_end - t_start) should be a multiple of dt_insert.")
    if t_end <= t_start:
        sys.exit("The t_end should be larger than the t_start.")
    if (abs(int(dt_insert/dt) - dt_insert/dt) >= dt) and (abs((int(dt_insert/dt)+1) - dt_insert/dt) >= dt):
        sys.exit("Variable dt_insert should be a multiple of time_step.")

    # Reading the inlet geometry and normal to the inlet condition
    coord_list = np.load(os.path.join(output_path, "inletPython.npy"))
    normal_inlet = np.load(os.path.join(output_path, "normalInletPython.npy"))

    # Initializing U and VOFw
    inlet_model = InletModel(config, coord_list, normal_inlet)

    probabilityShapes = [1.0]
    if np.sum(probabilityShapes) != 1.0:
        sys.exit('Vector "probabilityShapes" indicating the probability of occurrence of bubble shapes has not been defined correctly.')

    # The random generator selects:
    # - the bubble shape definition (shapeID)
    # - the center point of a bubble, both in
    # - inlet plane and in time
    # - the bubble mass
    n_time_insert = int((t_end-t_start)/dt_insert)
    logger.info(f"Between t_start {t_start} s and t_end {t_end} s, {n_time_insert} intervals of {dt_insert} s need to be defined.")
    for t_insert_idx in np.arange(n_time_insert):
        logger.info(f"Start bubble calculation for time interval {t_insert_idx}")
        mg_inserted = inlet_model.insert_bubbles(t_insert_idx)
        logger.info(f"\t Mass of inserted gas: {mg_inserted} kg. (Target mass: {mg_per_insert} kg).")
    logger.info(f"Inlet model iteration loop ended.")

    # Writing the profile to be used in OpenFOAM
    logger.info("Saving inlet profile to Python (numpy) npy-files.")
    np.save(os.path.join(output_path, "inletDefinition-U.npy"), inlet_model.velocity)
    np.save(os.path.join(output_path, "inletDefinition-VOFw.npy"), inlet_model.alpha)
    np.save(os.path.join(output_path, "inletDefinition-time.npy"), inlet_model.time)
    logger.info("Inlet profile saved in Python (numpy) npy-files.")

    # Check: convert to file compatible with ParaView to visualize the pre-inlet domain.
    logger.info("Saving inlet profile to CSV-files. ")
    files = [os.path.join(output_path, "inletDefinition-VOFw.csv"),
             os.path.join(output_path, "inletDefinition-Ux.csv"),
             os.path.join(output_path, "inletDefinition-Uy.csv"),
             os.path.join(output_path, "inletDefinition-Uz.csv"),
            ]
    toWrite = [inlet_model.alpha[:, :, 0], inlet_model.velocity[:, :, 0],
               inlet_model.velocity[:, :, 1], inlet_model.velocity[:, :, 2]]
    for fi in np.arange(len(files)):
        f = open(files[fi], 'w')
        f.write('x coord,y coord,z coord,value \n')
        for i in np.arange(len(coord_list[:, 0])):
            for j in np.arange(len(inlet_model.time)):
                coordPoint = np.array([coord_list[i, 1] - (velocity_bc * inlet_model.time[j]) * normal_inlet[0],
                                    coord_list[i, 2] - (velocity_bc * inlet_model.time[j]) * normal_inlet[1],
                                    coord_list[i, 3] - (velocity_bc * inlet_model.time[j]) * normal_inlet[2]])
                f.write(str(coordPoint[0]) + ',' + str(coordPoint[1]) + ',' + str(coordPoint[2]) + ',' + str(
                    toWrite[fi][i, j]) + '\n')
        f.close()
    logger.info("Inlet profile saved to CSV-files.")

    logger.info("========================End inlet_modelling========================")
