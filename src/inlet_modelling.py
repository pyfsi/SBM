from utils import np, sys, os, random, PI, SBM_OUTPUT, logging
logger = logging.getLogger(__name__)

class InletModel():
    def __init__(self, config, face_list, normal_inlet):

        # configuration
        self.t_start = float(config["cfd"]["start_time"])
        self.t_end = float(config["cfd"]["end_time"])
        self.timestep_size = float(config["cfd"]["delta_time"])
        self.block_size = float(config["sbm"]["t_unit"])
        self.density_gas = float(config["cfd"]["rho_g"])
        # self.rhol = float(config["cfd"]["rho_l"])
        self.mg_per_block = float(config["sbm"]["mass_gas"])
        self.tol_mg = float(config["sbm"]["mass_gas_tol"])
        self.velocity_bc = float(config["sbm"]["velocity"])
        self.mg_min = float(config["sbm"]["mass_bubble_min"])
        self.mg_max = float(config["sbm"]["mass_bubble_max"])
        self.intersect_boundary = config["sbm"]["intersect_boundary"]
        self.intersect_bubble = config["sbm"]["intersect_bubble"]

        # inlet geometry information
        self.face_list = face_list # ID - X - Y - Z - Area
        self.normal_inlet = normal_inlet

        # store previous bubble centers for next time step
        self.prev_centers = []

        # Initializing velocity and phase fraction arrays
        n_time_step = int((self.t_end - self.t_start)/self.timestep_size)+1
        self.velocity = np.ones([len(face_list), n_time_step, 3])
        self.velocity[:, :, 0] = self.velocity_bc * normal_inlet[0]
        self.velocity[:, :, 1] = self.velocity_bc * normal_inlet[1]
        self.velocity[:, :, 2] = self.velocity_bc * normal_inlet[2]
        self.alpha = np.ones([len(face_list), n_time_step, 1])
        self.time = np.linspace(self.t_start, self.t_end+self.timestep_size, n_time_step)

    def update(self, velocity, alpha):
        self.velocity = velocity
        self.alpha = alpha

    def sample_bubble_coord(self, face_idx_bounds, time_idx_bounds, bubble_mass_bounds):
        # sample cell index for spatial coordinates
        face_idx = random.randint(face_idx_bounds[0], face_idx_bounds[1])

        # sample cell index for temporal coordinates
        time_idx = random.randint(time_idx_bounds[0], time_idx_bounds[1])

        # sample bubble mass
        bubble_mass = random.uniform(bubble_mass_bounds[0], bubble_mass_bounds[1])

        return face_idx, time_idx, bubble_mass

    def define_bubble(self, face_idx, time_idx, block_idx, bubble_mass):

        t_start = self.t_start
        t_end = self.t_end
        timestep_size = self.timestep_size
        block_size = self.block_size
        density_gas = self.density_gas
        velocity_bc = self.velocity_bc
        intersect_boundary = self.intersect_boundary
        intersect_bubble = self.intersect_bubble
        face_list = self.face_list
        normal_inlet = self.normal_inlet
        time = self.time

        # temporary storage arrays (must be np.array)
        velocity_temp = np.array(self.velocity)
        alpha_temp = np.array(self.alpha)

        bubble_coord = face_list[face_idx, :]
        cell_time = time[time_idx]
        bubble_center = bubble_coord[1:4] - (velocity_bc * cell_time) * normal_inlet[:]

        # calculate gas radius assuming spherical bubble
        radius_gas = ((3.0*bubble_mass)/(4.0*PI*density_gas))**(1.0/3.0)

        rel_cell_time = cell_time - t_start - int((cell_time-t_start)/block_size) * block_size
        # Checks below prevents intersection with start and end of t_unit domain
        intersect_with_start = cell_time < radius_gas/velocity_bc
        intersect_with_end = cell_time > (t_end-radius_gas/velocity_bc)
        if intersect_with_start:
            return False, 0.0
        if intersect_with_end:
            return False, 0.0

        mg_defined = 0.0
        mg_bubble_wall = 0.0

        # get space and time index (i and j) of bounding box
        is_inside_bubble = np.linalg.norm(face_list[:, 1:4] - bubble_coord[1:4], axis=1) < radius_gas
        face_idx_in_bubble = is_inside_bubble.nonzero()[0]
        min_rel_time_idx_in_bubble = int((rel_cell_time - radius_gas/velocity_bc)/timestep_size)
        temp = (rel_cell_time + radius_gas/velocity_bc) // timestep_size
        max_rel_time_idx_in_bubble = int(temp) + 1

        # convert from relative time idx
        min_time_idx_in_bubble = block_idx * int(block_size / timestep_size) + min_rel_time_idx_in_bubble
        max_time_idx_in_bubble = block_idx * int(block_size / timestep_size) + max_rel_time_idx_in_bubble

        n_true_flag = 0
        for i in face_idx_in_bubble:
            cell_area = face_list[i, 4]
            for j in range(min_time_idx_in_bubble, max_time_idx_in_bubble):
                cell_coord = face_list[i, 1:4] - (velocity_bc * time[j]) * normal_inlet[:]
                if np.linalg.norm(cell_coord-bubble_center) < radius_gas:
                    if self.alpha[i, j, 0] == 1.0:
                        alpha_temp[i, j, 0] = 0.0
                        velocity_temp[i, j, :] = velocity_bc * normal_inlet[:]

                        mg_defined += cell_area * velocity_bc * timestep_size * density_gas
                        mg_bubble_wall += cell_area * velocity_bc * timestep_size * density_gas

                        n_true_flag += 1
                    elif intersect_bubble:
                        mg_bubble_wall += cell_area * velocity_bc * timestep_size * density_gas
                    else:
                        return False, 0.0

        if not(intersect_boundary):
            if mg_bubble_wall < (bubble_mass-np.average(face_list[:, 4])*velocity_bc*timestep_size):
                return False, 0.0

        # Save temporary data to class variable data
        self.update(velocity_temp, alpha_temp)

        return True, mg_defined

    def insert_bubbles(self, block_idx) -> float:
        # aliases
        mg_per_block = self.mg_per_block
        tol_mg = self.tol_mg
        face_list = self.face_list
        block_size = self.block_size
        timestep_size = self.timestep_size
        mg_min = self.mg_min
        mg_max = self.mg_max

        timesteps_per_block = int(block_size / timestep_size)
        min_time_at_blockidx = block_idx * timesteps_per_block
        max_time_at_blockidx = (block_idx + 1) * timesteps_per_block - 1
        n_inlet_faces = len(face_list) - 1

        iter = 0
        mg_inserted = 0.0
        # iterate until mass of inserted gas within tolerance of target mass
        while abs(mg_per_block-mg_inserted) > (tol_mg):
            # calculate bounds for bubble sampling
            face_idx_bounds = [0, n_inlet_faces]
            time_idx_bounds = [min_time_at_blockidx, max_time_at_blockidx]
            bubble_mass_bounds = [
                min((mg_min, mg_per_block - mg_inserted)),
                min((mg_max, mg_per_block - mg_inserted))
            ]

            face_idx, time_idx, bubble_mass = self.sample_bubble_coord(face_idx_bounds, time_idx_bounds, bubble_mass_bounds)
            is_bubble_defined, mg_defined = self.define_bubble(face_idx, time_idx, block_idx, bubble_mass)

            # print for debugging
            # residual = mg_per_block-mg_inserted
            # print(f"\t bubble_mass {bubble_mass:.3e} \t Residual {residual:.3e}")

            if is_bubble_defined:
                mg_inserted += mg_defined
                iter = 0
            else:
                iter = iter+1

            if iter > 1000:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations.")

        return mg_inserted

def inlet_modelling(config):
    logger.info("========================Start inlet_modelling========================")
    print(f"Running inlet_modelling")

    # configuration
    t_start = float(config["cfd"]["start_time"])
    t_end = float(config["cfd"]["end_time"])
    timestep_size = float(config["cfd"]["delta_time"])
    block_size = float(config["sbm"]["t_unit"])
    mg_per_block = float(config["sbm"]["mass_gas"])
    velocity_bc = float(config["sbm"]["velocity"])
    seed = config["sbm"].get("seed", None)
    output_path = os.path.join(os.getcwd(), SBM_OUTPUT)

    # set seed for random number generator
    if seed=="None":
        seed = None
    seed_msg = f"Using seed {seed} as {type(seed)} for random number generator"
    print(seed_msg)
    logging.info(seed_msg)
    random.seed(seed)

    if int((t_end-t_start)/block_size) != ((t_end-t_start)/block_size):
        raise RuntimeError("The desired time interval (t_end - t_start) should be a multiple of t_unit.")
    if t_end <= t_start:
        raise RuntimeError("The t_end should be larger than the t_start.")
    if (abs(int(block_size/timestep_size) - block_size/timestep_size) >= timestep_size) and (abs((int(block_size/timestep_size)+1) - block_size/timestep_size) >= timestep_size):
        raise RuntimeError("Variable block_size should be a multiple of time_step.")

    # Reading the inlet geometry and normal to the inlet condition
    face_list = np.load(os.path.join(output_path, "inletPython.npy"))
    normal_inlet = np.load(os.path.join(output_path, "normalInletPython.npy"))

    # Initialize inlet model class
    inlet_model = InletModel(config, face_list, normal_inlet)

    # The random generator selects:
    # - the center point of a bubble, both in inlet plane and in time
    # - the bubble mass
    n_blocks = int((t_end-t_start)/block_size) + 1
    logger.info(f"Between t_start {t_start} s and t_end {t_end} s, {n_blocks} intervals of {block_size} s need to be defined.")
    for block_idx in range(n_blocks - 1): # no insertion at last time step
        start_msg = f"Start bubble calculation for time interval {block_idx}"
        logger.info(start_msg)
        print(start_msg)
        mg_inserted = inlet_model.insert_bubbles(block_idx)
        logger.info(f"\t Mass of inserted gas: {mg_inserted} kg. (Target mass: {mg_per_block} kg).")
    logger.info(f"Inlet model iteration loop ended.")

    # Writing the profile to be used in OpenFOAM
    logger.info("Saving inlet profile to Python (numpy) npy-files.")
    np.save(os.path.join(output_path, "inletDefinition-U.npy"), inlet_model.velocity)
    np.save(os.path.join(output_path, "inletDefinition-VOFw.npy"), inlet_model.alpha)
    np.save(os.path.join(output_path, "inletDefinition-time.npy"), inlet_model.time)
    logger.info("Inlet profile saved in Python (numpy) npy-files.")

    # Check: convert to file compatible with ParaView to visualize the pre-inlet domain.
    logger.info("Saving inlet profile to CSV-files. ")
    csv_files = [
        os.path.join(output_path, "inletDefinition-VOFw.csv"),
        os.path.join(output_path, "inletDefinition-Ux.csv"),
        os.path.join(output_path, "inletDefinition-Uy.csv"),
        os.path.join(output_path, "inletDefinition-Uz.csv"),
    ]
    toWrite = [
        inlet_model.alpha[:, :, 0],
        inlet_model.velocity[:, :, 0],
        inlet_model.velocity[:, :, 1],
        inlet_model.velocity[:, :, 2]
    ]
    for fi, file_i in enumerate(csv_files):
        f = open(file_i, 'w')
        f.write('x coord,y coord,z coord,value \n')
        for i in np.arange(len(face_list[:, 0])):
            for j, time_j in enumerate(inlet_model.time):
                cell_coord = face_list[i, 1:4] - (velocity_bc * time_j) * normal_inlet[:]
                f.write(str(cell_coord[0]) + ',' + str(cell_coord[1]) + ',' + str(cell_coord[2]) + ',' +
                        str(toWrite[fi][i, j]) + '\n')
        f.close()
    logger.info("Inlet profile saved to CSV-files.")

    logger.info("========================End inlet_modelling========================")
