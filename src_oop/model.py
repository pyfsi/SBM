from utils import np, random, os
from utils import PI

class Model():
    def __init__(self, config, data, logger):
        # configuration
        self.density_gas = float(config["cfd"]["rho_g"])
        self.t_start = float(config["model"]["time"]["start"])
        self.t_end = float(config["model"]["time"]["end"])
        self.time_step = float(config["model"]["time"]["step"])
        self.block_size = float(config["model"]["time"]["block"])
        self.mg_per_block = float(config["model"]["mass_g"]["per_block"])
        self.mg_tol = float(config["model"]["mass_g"]["tol"])
        self.mg_min = float(config["model"]["mass_g"]["min"])
        self.mg_max = float(config["model"]["mass_g"]["max"])
        self.velocity_bc = float(config["model"]["velocity"])
        self.intersect_boundary = config["model"]["intersect_boundary"]
        self.intersect_bubble = config["model"]["intersect_bubble"]
        self.seed = str(config["model"].get("seed", None))
        self.output_path = str(config.get("_output_path"))

        # reference to data storage
        self.data = data

        # logger
        self.logger = logger

    def initialize(self):
        '''initialize data storage for physical variables at the inlet'''
        self.time = np.arange(self.t_start, self.t_end+self.time_step, self.time_step)

        self.n_timesteps = len(self.time)
        self.n_faces = len(self.data["inlet_faces"])
        self.alpha = np.ones([self.n_faces, self.n_timesteps, 1], dtype=np.float64)
        self.velocity = np.ones([self.n_faces, self.n_timesteps, 3], dtype=np.float64)
        self.velocity[:, :, 0] = self.velocity_bc * self.data["inlet_normal"][0]
        self.velocity[:, :, 1] = self.velocity_bc * self.data["inlet_normal"][1]
        self.velocity[:, :, 2] = self.velocity_bc * self.data["inlet_normal"][2]

    def run(self):
        # read configuration
        t_start = self.t_start
        t_end = self.t_end
        block_size = self.block_size
        mg_per_block = self.mg_per_block

        # inlet modelling
        self._set_seed()
        n_blocks = int((t_end-t_start)/block_size)
        self.logger.info(f"Discretized {n_blocks} intervals of {block_size} s between t_start {t_start} s and t_end {t_end} s.")
        for block_idx in range(n_blocks):
            start_msg = f"Start bubble calculation for time interval {block_idx}"
            self.logger.info(start_msg)
            print(start_msg)
            mass_inserted = self._insert_bubbles(block_idx)
            self.logger.info(f"\t Mass of inserted gas: {mass_inserted} kg. (Target mass: {mg_per_block} kg).")
        self.logger.info(f"Inlet model iteration loop ended.")
        self._pass_data()

        # save
        self._save_inlet_data()
        self._save_csv()

    def get_data(self):
        return self.time, self.alpha, self.velocity

    # ===== Protected functions =====
    def _convert_relative_time_idx(self, block_idx: int, val: int) -> int:
            '''
            Convert time index relative to time blocks to absolute time index.
            Args:
                val: relative time index

            Returns:
                Returns the absolute time index
            '''
            return block_idx * int(self.block_size / self.time_step) + val

    def _set_seed(self):
        # set seed for random number generator
        seed = self.seed
        if seed.lower()=="none":
            seed = None
        seed_msg = f"Using seed {seed} for random number generator"
        print(seed_msg)
        random.seed(seed)

    def _sample_random(self, fid_bounds: list, tid_bounds: list, mass_bounds:list):
        '''
        Sample bubble parameters from predefined statistical distributions.
        Args:
            fid_bounds: lower and upper bounds for face index
            tid_bounds: lower and upper bounds for time index
            mass_bounds: lower and upper bounds for mass of one bubble

        Returns:
            Returns face_idx, time_idx, mass
        '''
        # sample cell index for spatial coordinates
        face_idx = random.randint(fid_bounds[0], fid_bounds[1])

        # sample time index for temporal coordinates
        time_idx = random.randint(tid_bounds[0], tid_bounds[1])

        # sample bubble mass
        mass_sample = random.uniform(mass_bounds[0], mass_bounds[1])

        return face_idx, time_idx, mass_sample

    def _define_bubble(self, face_idx, time_idx, mass_sample, block_idx):
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
        time_step = self.time_step
        block_size = self.block_size
        density_gas = self.density_gas
        velocity_bc = self.velocity_bc
        intersect_boundary = self.intersect_boundary
        intersect_bubble = self.intersect_bubble
        face_list = self.data["inlet_faces"]
        normal_inlet = self.data["inlet_normal"]
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
        min_rel_time_idx_in_radius = int((rel_cell_time - radius_gas/velocity_bc) / time_step)
        temp = (rel_cell_time + radius_gas/velocity_bc) // time_step
        max_rel_time_idx_in_radius = int(temp) + 1

        # convert from relative time idx
        min_time_idx_in_radius = self._convert_relative_time_idx(block_idx, min_rel_time_idx_in_radius)
        max_time_idx_in_radius = self._convert_relative_time_idx(block_idx, max_rel_time_idx_in_radius)
        time_idx_in_radius = np.arange(min_time_idx_in_radius, max_time_idx_in_radius)

        # create 3d tensor to describe the relative cell positions w.r.t. velocity times time
        face_list_extended = np.array([face_list[face_idx_in_radius, 1:4]] * len(time_idx_in_radius))
        time_velocity_product = np.tensordot(time[time_idx_in_radius], velocity_bc * normal_inlet[:], axes=0)
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
        avg_gas_mass_per_cell = np.average(face_list[:, 4]) * velocity_bc * time_step * density_gas
        cell_area_inside_bubble = np.sum(face_list[face_idx_in_bubble, 4])
        bubble_mass_defined = cell_area_inside_bubble * velocity_bc * time_step * density_gas

        # return False if defined bubble mass is smaller than expected
        # and bubble is not allowed to intersect boundary
        if not intersect_boundary:
            if (mass_sample-bubble_mass_defined) > avg_gas_mass_per_cell:
                return False, 0.0

        # set alpha and velocity fields and defined gas mass
        self.alpha[face_idx_in_bubble, time_idx_in_bubble, 0] = 0.0
        self.velocity[face_idx_in_bubble, time_idx_in_bubble, :] = velocity_bc * normal_inlet[:]
        mass_bubble_cells = cell_area_inside_bubble * density_gas * velocity_bc * time_step

        return True, mass_bubble_cells

    def _insert_bubbles(self, block_idx: int) -> float:
        '''
        Iteratively insert bubble to the block until target mass has been reached.

        Args:
            block_idx: index for insertion block

        Returns:
            mass_inserted: total mass of cells enclosed by defined bubbles.
        '''
        # aliases
        mass_per_block = self.mg_per_block
        mg_tol = self.mg_tol
        face_list = self.data["inlet_faces"]
        block_size = self.block_size
        time_step = self.time_step
        mass_lower_bound = self.mg_min
        mass_upper_bound = self.mg_max

        # calc time indices
        timesteps_per_block = int(block_size / time_step)
        min_time_at_blockidx = self._convert_relative_time_idx(block_idx, 0)
        max_time_at_blockidx = self._convert_relative_time_idx(block_idx, timesteps_per_block - 1)

        iter = 0
        mass_inserted = 0.0
        # iterate until mass of inserted gas is within tolerance of target mass
        while abs(mass_per_block-mass_inserted) > (mg_tol):
            # calculate bounds of sample space for bubble parameters
            face_idx_bounds = [0, len(face_list) - 1]
            time_idx_bounds = [min_time_at_blockidx, max_time_at_blockidx]
            mass_bounds = [
                min((mass_lower_bound, mass_per_block - mass_inserted)),
                min((mass_upper_bound, mass_per_block - mass_inserted))
            ]

            face_idx, time_idx, mass_sample = self._sample_random(face_idx_bounds, time_idx_bounds, mass_bounds)
            is_bubble_defined, mass_bubble_cells = self._define_bubble(face_idx, time_idx, mass_sample, block_idx)

            # if self.plotter:
            #     self.plotter.update_sample_distribution(mass_sample)

            if is_bubble_defined:
                mass_inserted += mass_bubble_cells
                iter = 0

                # # run plotter
                # if self.plotter:
                #     self.plotter.update_data(face_idx, time_idx)
                #     self.plotter.update_distribution(mass_bubble_cells)
                #     self.plotter.update_residual(mass_per_block, mass_inserted)
                #     self.plotter.plot_bubble_insertion()

                # log face_idx, time_idx, and bubble mass
                self.logger.info(f"\t\t inserted at face_idx={face_idx}, time_idx={time_idx}, m_b={mass_sample}")
            else:
                iter = iter+1

            if iter > 1000:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations.")

        return mass_inserted

    def _pass_data(self):
        self.data["time"] = self.time
        self.data["alpha"] = self.alpha
        self.data["velocity"] = self.velocity

    def _save_inlet_data(self):
        self.logger.info("Saving inlet profile to npy-files.")
        np.save(os.path.join(self.output_path, "inlet_velocity.npy"), self.velocity)
        np.save(os.path.join(self.output_path, "inlet_alpha.npy"), self.alpha)
        np.save(os.path.join(self.output_path, "inlet_time.npy"), self.time)
        self.logger.info("Inlet profile saved in npy-files.")

    def _save_csv(self):
        # print csv to visualize pre-inlet domain
        self.logger.info("Saving inlet profile to csv-files.")
        csv_file_path = os.path.join(self.output_path, "inlet_data.csv")
        inlet_all_variable = np.concatenate((self.alpha[:,:,:], self.velocity[:,:,:]), axis=2)
        csv_header = "x_coord,y_coord,z_coord,alpha,velocity_x,velocity_y,velocity_z"

        # get cell coordinates in x,y,z space
        n_timesteps = self.n_timesteps
        n_faces = self.n_faces
        face_list_extended = np.array([self.data["inlet_faces"][:, 1:4]] * n_timesteps)
        time_velocity_product = np.tensordot(self.time[:], self.velocity_bc * self.data["inlet_normal"][:], axes=0)
        cell_coords = face_list_extended[:, :, :] - time_velocity_product[:, None, : ]
        cell_coords = np.swapaxes(cell_coords, 0, 1)
        cell_coords = np.reshape(cell_coords, (n_timesteps * n_faces, -1), order='C')

        # save csv
        inlet_var_reshaped = np.reshape(inlet_all_variable, (n_timesteps * n_faces, -1), order="C")
        inlet_ds = np.concatenate((cell_coords, inlet_var_reshaped), axis=1)
        np.savetxt(csv_file_path, inlet_ds, fmt='%.6e',
                    header=csv_header, delimiter=",", comments='')
        self.logger.info("Inlet profile saved to csv-files.")