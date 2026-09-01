from utils import np, sys, os, random, logging
from utils import PI, SBM_OUTPUT
from .plotter import Plotter

logger = logging.getLogger(__name__)

class InletModel():
    def __init__(self, config, face_list, normal_inlet, plotter):

        # configuration
        self.t_start = float(config["cfd"]["start_time"])
        self.t_end = float(config["cfd"]["end_time"])
        self.timestep_size = float(config["cfd"]["delta_time"])
        self.block_size = float(config["sbm"]["t_unit"])
        self.density_gas = float(config["cfd"]["rho_g"])
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

        # initializing velocity and phase fraction arrays
        self.time = np.arange(self.t_start, self.t_end+self.timestep_size, self.timestep_size)
        n_time_step = len(self.time)
        self.velocity = np.ones([len(face_list), n_time_step, 3], dtype=np.float64)
        self.velocity[:, :, 0] = self.velocity_bc * normal_inlet[0]
        self.velocity[:, :, 1] = self.velocity_bc * normal_inlet[1]
        self.velocity[:, :, 2] = self.velocity_bc * normal_inlet[2]
        self.alpha = np.ones([len(face_list), n_time_step, 1], dtype=np.float64)

        # plotter
        self.plotter = plotter

    def _convert_relative_time_idx(self, block_idx: int, val: int) -> int:
        '''
        Convert time index relative to time blocks to absolute time index.
        Args:
            val: relative time index

        Returns:
            Returns the absolute time index
        '''
        return block_idx * int(self.block_size / self.timestep_size) + val

    def sample_random(self, face_idx_bounds: list, time_idx_bounds: list, mass_bounds: list):
        '''
        Sample bubble parameters from predefined statistical distributions.
        Args:
            face_idx_bounds: lower and upper bounds for face index
            time_idx_bounds: lower and upper bounds for time index
            mass_bounds: lower and upper bounds for mass of one bubble

        Returns:
            Returns face_idx, time_idx, mass

        '''
        # sample cell index for spatial coordinates
        face_idx = random.randint(face_idx_bounds[0], face_idx_bounds[1])

        # sample time index for temporal coordinates
        time_idx = random.randint(time_idx_bounds[0], time_idx_bounds[1])

        # sample bubble mass
        mass_sample = random.uniform(mass_bounds[0], mass_bounds[1])

        return face_idx, time_idx, mass_sample

    def define_bubble(self, face_idx, time_idx, mass_sample, block_idx):
        '''
        Set volume of fluid fraction and velocity cells within a bubble to their prescribed values.

        Args:
            face_idx: index for inlet faces array
            time_idx: index for time array
            mass_sample: mass of one bubble
            block_idx: index for insertion block

        Returns:
            is_bubble_defined: boolean describing validity of bubble
            mass_bubble_cells: mass of cells enclosed by the bubble
        '''

        # alias
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

        bubble_coord = face_list[face_idx, :] # ID - X - Y - Z - area
        bubble_time = time[time_idx]
        bubble_center = bubble_coord[1:4] - (velocity_bc * bubble_time) * normal_inlet[:] # X - Y - Z

        # calculate gas radius assuming spherical bubble
        radius_gas = ((3.0*mass_sample)/(4.0*PI*density_gas))**(1.0/3.0)

        rel_cell_time = bubble_time - t_start - int((bubble_time-t_start)/block_size) * block_size
        # Checks below prevents intersection with start and end of t_unit domain
        intersect_with_start = bubble_time < radius_gas/velocity_bc
        intersect_with_end = bubble_time > (t_end-radius_gas/velocity_bc)
        if intersect_with_start:
            return False, 0.0
        if intersect_with_end:
            return False, 0.0

        # get space and time index (i and j) of bounding box
        is_inside_radius = np.linalg.norm(face_list[:, 1:4] - bubble_coord[1:4], axis=1) < radius_gas
        face_idx_in_radius = is_inside_radius.nonzero()[0]
        min_rel_time_idx_in_radius = int((rel_cell_time - radius_gas/velocity_bc) / timestep_size)
        temp = (rel_cell_time + radius_gas/velocity_bc) // timestep_size
        max_rel_time_idx_in_radius = int(temp) + 1

        # convert from relative time idx
        min_time_idx_in_radius = self._convert_relative_time_idx(block_idx, min_rel_time_idx_in_radius)
        max_time_idx_in_radius = self._convert_relative_time_idx(block_idx, max_rel_time_idx_in_radius)
        time_idx_in_radius = np.arange(min_time_idx_in_radius, max_time_idx_in_radius)

        # create 3d tensor to describe the relative cell positions w.r.t. velocity times time
        face_list_extended = np.array([face_list[face_idx_in_radius, 1:4]] * len(time_idx_in_radius))
        time_velocity_product = np.tensordot(time[time_idx_in_radius], velocity_bc*normal_inlet[:], axes=0)
        cell_coords = face_list_extended[:,:,:] - time_velocity_product[:,None,:]
        cell_coords = np.swapaxes(cell_coords, 0, 1) # [faces, timesteps, xyz]

        # get boolean field if cell is inside bubble
        displacement = cell_coords[:,:,:]-bubble_center[None,None,:]
        distance_sqr = np.sum(displacement*displacement, axis=2)
        is_cell_inside_bubble = distance_sqr < radius_gas * radius_gas
        face_idx_in_bubble = face_idx_in_radius[is_cell_inside_bubble.nonzero()[0]]
        time_idx_in_bubble = time_idx_in_radius[is_cell_inside_bubble.nonzero()[1]]

        # return False if alpha inside to-be-defined bubble is already zero
        alpha_inside_bubble = self.alpha[face_idx_in_bubble, time_idx_in_bubble, 0]
        if not intersect_bubble:
            if np.any(alpha_inside_bubble==0.0):
                return False, 0.0

        # calculate defined bubble mass
        avg_gas_mass_per_cell = np.average(face_list[:, 4]) * velocity_bc * timestep_size * density_gas
        cell_area_inside_bubble = np.sum(face_list[face_idx_in_bubble, 4])
        bubble_mass_defined = cell_area_inside_bubble * velocity_bc * timestep_size * density_gas

        # return False if defined bubble mass is smaller than expected
        # and bubble is not allowed to intersect boundary
        if not intersect_boundary:
            if (mass_sample-bubble_mass_defined) > avg_gas_mass_per_cell:
                return False, 0.0

        # set alpha and velocity fields and defined gas mass
        self.alpha[face_idx_in_bubble, time_idx_in_bubble, 0] = 0.0
        self.velocity[face_idx_in_bubble, time_idx_in_bubble, :] = velocity_bc * normal_inlet[:]
        mass_bubble_cells = cell_area_inside_bubble * density_gas * velocity_bc * timestep_size

        return True, mass_bubble_cells

    def insert_bubbles(self, block_idx: int) -> float:
        '''
        Iteratively insert bubble to the block until target mass has been reached.

        Args:
            block_idx: index for insertion block

        Returns:
            mass_inserted: total mass of cells enclosed by defined bubbles.
        '''
        # aliases
        mass_per_block = self.mg_per_block
        tol_mg = self.tol_mg
        face_list = self.face_list
        block_size = self.block_size
        timestep_size = self.timestep_size
        mass_lower_bound = self.mg_min
        mass_upper_bound = self.mg_max

        # calc time indices
        timesteps_per_block = int(block_size / timestep_size)
        min_time_at_blockidx = self._convert_relative_time_idx(block_idx, 0)
        max_time_at_blockidx = self._convert_relative_time_idx(block_idx, timesteps_per_block - 1)

        iter = 0
        mass_inserted = 0.0
        # iterate until mass of inserted gas is within tolerance of target mass
        while abs(mass_per_block-mass_inserted) > (tol_mg):
            # calculate bounds of sample space for bubble parameters
            face_idx_bounds = [0, len(face_list) - 1]
            time_idx_bounds = [min_time_at_blockidx, max_time_at_blockidx]
            mass_bounds = [
                min((mass_lower_bound, mass_per_block - mass_inserted)),
                min((mass_upper_bound, mass_per_block - mass_inserted))
            ]

            face_idx, time_idx, mass_sample = self.sample_random(face_idx_bounds, time_idx_bounds, mass_bounds)
            is_bubble_defined, mass_bubble_cells = self.define_bubble(face_idx, time_idx, mass_sample, block_idx)

            if self.plotter:
                self.plotter.update_sample_distribution(mass_sample)

            if is_bubble_defined:
                mass_inserted += mass_bubble_cells
                iter = 0

                # run plotter
                if self.plotter:
                    self.plotter.update_data(face_idx, time_idx)
                    self.plotter.update_distribution(mass_bubble_cells)
                    self.plotter.update_residual(mass_per_block, mass_inserted)
                    self.plotter.plot_bubble_insertion()

                # log face_idx, time_idx, and bubble mass
                logger.info(f"\t\t inserted at face_idx={face_idx}, time_idx={time_idx}, m_b={mass_sample}")
            else:
                iter = iter+1

            if iter > 1000:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations.")

        return mass_inserted

def inlet_modelling(config):
    print(f"Running inlet_modelling")

    # read configuration
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

    # reading inlet geometry and normal to the inlet condition
    face_list = np.load(os.path.join(output_path, "inletPython.npy"))
    normal_inlet = np.load(os.path.join(output_path, "normalInletPython.npy"))

    # initialize distribution plotter
    activate_plotter = config["activate_plotter"]
    plotter = None
    if activate_plotter:
        plotter = Plotter(config)

    # initialize inlet model class
    inlet_model = InletModel(config, face_list, normal_inlet, plotter)

    # iterate through every insertion blocks
    n_blocks = int((t_end-t_start)/block_size)
    logger.info(f"Discretized {n_blocks} intervals of {block_size} s between t_start {t_start} s and t_end {t_end} s.")
    for block_idx in range(n_blocks):
        start_msg = f"Start bubble calculation for time interval {block_idx}"
        logger.info(start_msg)
        print(start_msg)
        mass_inserted = inlet_model.insert_bubbles(block_idx)
        logger.info(f"\t Mass of inserted gas: {mass_inserted} kg. (Target mass: {mg_per_block} kg).")
    logger.info(f"Inlet model iteration loop ended.")

    # save plot
    if activate_plotter:
        plotter.save_plot(output_path)

    # Write alpha and velocity profiles
    logger.info("Saving inlet profile to Python (numpy) npy-files.")
    np.save(os.path.join(output_path, "inletDefinition-U.npy"), inlet_model.velocity)
    np.save(os.path.join(output_path, "inletDefinition-VOFw.npy"), inlet_model.alpha)
    np.save(os.path.join(output_path, "inletDefinition-time.npy"), inlet_model.time)
    logger.info("Inlet profile saved in Python (numpy) npy-files.")

    # print csv to visualize pre-inlet domain
    logger.info("Saving inlet profile to CSV-files. ")
    csv_file_path = os.path.join(output_path, "inletModelData.csv")
    inlet_all_variable = np.concatenate((inlet_model.alpha[:,:,:], inlet_model.velocity[:,:,:]), axis=2)
    csv_header = "x_coord,y_coord,z_coord,alpha,velocity_x,velocity_y,velocity_z"

    # get cell coordinates in x,y,z space
    n_timesteps = len(inlet_model.time)
    n_faces = len(face_list)
    face_list_extended = np.array([face_list[:, 1:4]] * n_timesteps)
    time_velocity_product = np.tensordot(inlet_model.time[:], velocity_bc*normal_inlet[:], axes=0)
    cell_coords = face_list_extended[:, :, :] - time_velocity_product[:, None, : ]
    cell_coords = np.swapaxes(cell_coords, 0, 1)
    cell_coords = np.reshape(cell_coords, (n_timesteps*n_faces, -1), order='C')

    # save csv
    inlet_var_reshaped = np.reshape(inlet_all_variable, (n_timesteps*n_faces, -1), order="C")
    inlet_ds = np.concatenate((cell_coords, inlet_var_reshaped), axis=1)
    np.savetxt(csv_file_path, inlet_ds, fmt='%.6e',
               header=csv_header, delimiter=",", comments='')
    logger.info("Inlet profile saved to CSV-files.")
