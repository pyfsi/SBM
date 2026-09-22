from utils import np, os
from utils import PI
from .sample_generator import SampleGenerator

class Model():
    MAX_INSERT_ITER = 1000

    def __init__(self, config, data, logger):
        # configuration
        self.density_gas = float(config["cfd"]["rho_g"])
        self.time_start = float(config["model"]["time"]["start"])
        self.time_end = float(config["model"]["time"]["end"])
        self.time_step = float(config["model"]["time"]["step"])
        self.time_block = float(config["model"]["time"]["block"])
        self.mg_per_block = float(config["model"]["mass_g"]["per_block"])
        self.mg_tol = float(config["model"]["mass_g"]["tol"])
        self.mg_min = float(config["model"]["mass_g"]["min"])
        self.mg_max = float(config["model"]["mass_g"]["max"])
        self.velocity_bc = float(config["model"]["velocity"])
        self.intersect_boundary = config["model"]["intersect_boundary"]
        self.intersect_bubble = config["model"]["intersect_bubble"]
        self.seed = str(config["model"].get("seed", None))

        self.output_path = str(config.get("_output_path"))
        self.timesteps_per_block = int(config.get("_timesteps_per_block"))
        self.buffer_size = int(config.get("_buffer_size"))

        # reference to data storage
        self.inlet_data = data
        self.inlet_data.faces = self.inlet_data.faces
        self.inlet_normal = self.inlet_data.normal

        # init generator for bubble definition
        self.generator = SampleGenerator(config,)

        # logger
        self.logger = logger

    def initialize(self, block_idx):
        # absolute time index and time array
        abs_time_idx_block_start = self.timesteps_per_block * block_idx
        abs_time_idx_block_end = self.timesteps_per_block * (block_idx+1) + self.buffer_size
        self.abs_time_idx = np.arange(abs_time_idx_block_start, abs_time_idx_block_end, 1)
        self.time = self.abs_time_idx * self.time_step + self.time_start
        self.n_timesteps = len(self.time)

        self.alpha = np.ones([self.n_faces, self.n_timesteps, 1], dtype=np.float64)
        self.velocity = np.ones([self.n_faces, self.n_timesteps, 3], dtype=np.float64)
        self.velocity[:,:,:] = self.velocity_bc * self.inlet_data.normal[None, None, :]

        # print log
        start_msg = f"Start bubble calculation for time interval {block_idx}"
        self.logger.info(start_msg)
        print(start_msg)

    def initialize_buffer(self):
        self.n_faces = len(self.inlet_data.faces)

        # buffer to pass alpha and velocity to next block iteration
        self.buffer = np.ones([self.n_faces, self.buffer_size, 4], dtype=np.float64)
        self.buffer[:,:,1:4] = self.velocity_bc * self.inlet_data.normal[None, None, :]

    def run(self):
        self._read_buffer()

        mass_inserted = self._insert_bubbles()
        inserted_mass_msg = f"\t Mass of inserted gas: {mass_inserted} kg. (Target mass: {self.mg_per_block} kg)."
        self.logger.info(inserted_mass_msg)

        self.inlet_data.store_block_data(self.alpha, self.velocity)

        self._store_buffer()

    # ===== Protected functions =====

    def _define_bubble(self, sample):
        '''
        Set volume of fluid fraction and velocity cells within a bubble to their prescribed values.

        Args:
            sample:
            block_idx: index for insertion block

        Returns:
            is_bubble_defined: boolean describing validity of bubble
            mass_bubble_cells: mass of cells enclosed by the bubble
        '''

        # alias
        time_start = self.time_start
        time_end = self.time_end
        time_step = self.time_step
        time_block = self.time_block
        density_gas = self.density_gas
        velocity_bc = self.velocity_bc
        intersect_boundary = self.intersect_boundary
        intersect_bubble = self.intersect_bubble
        # arrays
        time = self.time
        face_list = self.inlet_data.faces
        normal_inlet = self.inlet_data.normal

        # sample
        face_idx = sample["face"]
        time_idx = sample["time"]
        mass_sample = sample["mass"]

        bubble_coord = face_list[face_idx, :] # ID - X - Y - Z - area
        bubble_time = time[time_idx]
        bubble_center = bubble_coord[1:4] - (velocity_bc * bubble_time) * normal_inlet[:] # X - Y - Z

        # calculate gas radius assuming spherical bubble
        radius_bubble = ((3.0*mass_sample)/(4.0*PI*density_gas))**(1.0/3.0)
        rel_cell_time = bubble_time - time_start - int((bubble_time-time_start)/time_block) * time_block

        # Checks below prevents intersection with start and end of the insertion block
        intersect_with_start = bubble_time < (time_start+radius_bubble/velocity_bc)
        intersect_with_end = bubble_time > (time_end-radius_bubble/velocity_bc)
        if intersect_with_start:
            return False, 0.0
        if intersect_with_end:
            return False, 0.0

        # get face and time index (i and j) of bounding box
        is_inside_radius = np.linalg.norm(face_list[:, 1:4] - bubble_coord[1:4], axis=1) < radius_bubble
        face_idx_in_radius = is_inside_radius.nonzero()[0]
        min_rel_time_idx_in_radius = int((rel_cell_time - radius_bubble/velocity_bc) / time_step)
        temp = (rel_cell_time + radius_bubble/velocity_bc) // time_step
        max_rel_time_idx_in_radius = int(temp) + 1

        # convert from relative time idx
        time_idx_in_radius = np.arange(min_rel_time_idx_in_radius, max_rel_time_idx_in_radius) # relative time idx

        # create 3d tensor to describe the relative cell positions w.r.t. velocity times time
        face_list_extended = np.array([face_list[face_idx_in_radius, 1:4]] * len(time_idx_in_radius))
        time_velocity_product = np.tensordot(time[time_idx_in_radius], velocity_bc * normal_inlet[:], axes=0)
        cell_coords = face_list_extended[:,:,:] - time_velocity_product[:,None,:]
        cell_coords = np.swapaxes(cell_coords, 0, 1) # [faces, timesteps, xyz]

        # get boolean field if cell is inside bubble
        displacement = cell_coords[:,:,:]-bubble_center[None,None,:]
        distance_sqr = np.sum(displacement*displacement, axis=2)
        is_cell_inside_bubble = distance_sqr < radius_bubble * radius_bubble
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

    def _insert_bubbles(self) -> float:
        '''
        Iteratively insert bubble to the block until target mass has been reached.

        Args:
            -

        Returns:
            mass_inserted: total mass of cells enclosed by defined bubbles.
        '''
        # aliases
        mass_per_block = self.mg_per_block
        mg_tol = self.mg_tol
        face_list = self.inlet_data.faces
        mass_lower_bound = self.mg_min
        mass_upper_bound = self.mg_max

        # calc time indices
        timesteps_per_block = self.timesteps_per_block
        min_timeidx = 0
        max_timeidx = timesteps_per_block - 1

        iter = 0
        mass_inserted = 0.0
        # iterate until mass of inserted gas is within tolerance of target mass
        while abs(mass_per_block-mass_inserted) > (mg_tol):
            # calculate bounds of sample space for bubble parameters
            face_idx_bounds = [0, len(face_list) - 1]
            time_idx_bounds = [min_timeidx, max_timeidx]
            mass_bounds = [
                min((mass_lower_bound, mass_per_block - mass_inserted)),
                min((mass_upper_bound, mass_per_block - mass_inserted))
            ]

            sample = self.generator.sample_random(face_idx_bounds, time_idx_bounds, mass_bounds)
            is_bubble_defined, mass_bubble_cells = self._define_bubble(sample)

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
                self.logger.info(f"\t\t inserted at \
                                 face_idx={sample["face"]}, \
                                 time_idx={sample["time"]}, \
                                 m_b={sample["mass"]}")
            else:
                iter = iter+1

            if iter > self.MAX_INSERT_ITER:
                raise RuntimeError("inlet_modelling took longer than 1000 iterations.")

        return mass_inserted

    # === conversion ===
    def _convert_to_time_idx(self, block_idx: int, rel_time_idx):
        '''
        Convert time index relative to time blocks to absolute time index.
        Args:
            rel_time_idx: relative time index

        Returns:
            the absolute time index
        '''
        return block_idx * self.timesteps_per_block + rel_time_idx

    def _convert_to_rel_time_idx(self, block_idx: int, time_idx):
        '''
        Convert time index relative to time blocks to absolute time index.
        Args:
            time_idx: relative time index

        Returns:
            the absolute time index
        '''
        return time_idx - block_idx * self.timesteps_per_block

    # === buffer operations ===
    def _read_buffer(self):
        b_size = self.buffer_size
        self.alpha[:, :b_size, 0] = self.buffer[:,:,0]
        self.velocity[:, :b_size, :] = self.buffer[:,:,1:4]

    def _store_buffer(self):
        b_size = self.buffer_size
        self.buffer[:,:,0] = self.alpha[:, -b_size:, 0]
        self.buffer[:,:,1:4] = self.velocity[:, -b_size:, :]

    # === save ===
    def _save_npy(self):
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
        face_list_extended = np.array([self.inlet_data.faces[:, 1:4]] * n_timesteps)
        time_velocity_product = np.tensordot(self.time[:], self.velocity_bc * self.inlet_data.normal[:], axes=0)
        cell_coords = face_list_extended[:, :, :] - time_velocity_product[:, None, : ]
        cell_coords = np.swapaxes(cell_coords, 0, 1)
        cell_coords = np.reshape(cell_coords, (n_timesteps * n_faces, -1), order='C')

        # save csv
        inlet_var_reshaped = np.reshape(inlet_all_variable, (n_timesteps * n_faces, -1), order="C")
        inlet_ds = np.concatenate((cell_coords, inlet_var_reshaped), axis=1)
        np.savetxt(csv_file_path, inlet_ds, fmt='%.6e',
                    header=csv_header, delimiter=",", comments='')
        self.logger.info("Inlet profile saved to csv-files.")