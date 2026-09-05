from utils import os, shutil, re, np
from utils import modulo

# SBM components
from .logger import Logger
from .preprocessor import Preprocessor
from .reader import Reader
from .model import Model
from .writer import Writer

class SBM():
    """Main class of the Synthetic Bubble Model."""
    def __init__(self, config):
        # config dictionary
        self.config = config
        self.time_start = float(config["model"]["time"]["start"])
        self.time_end = float(config["model"]["time"]["end"])
        self.time_step = float(config["model"]["time"]["step"])
        self.time_block = float(config["model"]["time"]["block"])
        self.velocity_bc = float(config["model"]["velocity"])
        self.inlet_name = str(config["cfd"]["inlet_name"])

        # data storage per insertion block
        self.data = {}
        self.data["inlet_faces"] = None
        self.data["inlet_normal"] = None
        self.data["time"] = None
        self.data["alpha"] = None
        self.data["velocity"] = None

        # create 'sbm_files'
        self._make_output_dir()
        self._add_time_settings()

    def check_case(self):
        # check SBM config
        self._check_time_settings()
        self._check_mass_settings()
        self._check_modules()

    def purge_boundary_data(self):
        purge_boundary_data = self.config["settings"]["purge_boundary_data"]
        self.config["_boundary_data_path"] = os.path.join(self.cwd, "constant", "boundaryData")

        if purge_boundary_data:
            if os.path.exists(self.config["_boundary_data_path"]):
                shutil.rmtree(self.config["_boundary_data_path"])
        else:
            raise RuntimeError("'boundaryData' directory already exists.")

    def initialize(self):
        config = self.config
        data = self.data

        # === initialize sbm components ===
        # diagnostics
        # self.profiler = Profiler(config)
        self.logger = Logger(config)

        # preprocessor
        self.preprocessor = Preprocessor(config)

        # reader
        self.reader = Reader(config, data, self.logger)

        # inlet modelling
        self.model = Model(config, data, self.logger)

        # writer
        self.writer = Writer(config, data, self.logger)
        self.writer.check()

        # plotter
        self.plotter = None

    def run(self):
        with self.logger.function_call(name="preprocessor"):
            self.preprocessor.initialize()
            self.preprocessor.run()

        with self.logger.function_call(name="reader"):
            self.reader.initialize()
            self.reader.run()

        self._initialize_cells()
        self._prepare_boundary_data_dir()
        self._iterate()

        # save
        self._save_csv()

    def finalize(self):
        # attributes
        cwd = self.cwd
        output_path = self.config["_output_path"]
        openfoam_type = self.config["_openfoam_type"]

        if openfoam_type=="org":
            filenames = ["Ccx", "Ccy", "Ccz", "area"]
        elif openfoam_type=="com":
            filenames = ["Cx", "Cy", "Cz", "area"]
        else:
            filenames = []

        # move files to 'sbm_files'
        for file in filenames:
            target_area_path = os.path.join(output_path, file)
            if os.path.exists(target_area_path):
                os.remove(target_area_path)
            area_path = os.path.join(cwd, "0", file)
            if os.path.exists(area_path):
                shutil.move(area_path, output_path)

    # == protected functions ==
    def _make_output_dir(self):
        # make sbm output dir
        self.cwd = os.getcwd()
        sbm_output_path = os.path.join(self.cwd, "sbm_files")
        # remove previous sbm output files
        if os.path.exists(sbm_output_path):
            shutil.rmtree(sbm_output_path)
            print(f"Deleting previous SBM output files.")
        os.mkdir(sbm_output_path)

        # === add attributes to config ===
        self.config["_output_path"] = sbm_output_path
        is_cfd_solver_openfoam = (self.config["packages"]["cfd_program"].lower() == "openfoam")
        if is_cfd_solver_openfoam:
            self.config["_openfoam_type"] = self._get_openfoam_type()

    def _add_time_settings(self):
        time_start = self.time_start
        time_end = self.time_end
        time_step = self.time_step
        time_block = self.time_block

        timesteps_per_block = int(time_block/time_step)
        self.config["_timesteps_per_block"] = timesteps_per_block
        buffer_size = int(timesteps_per_block*0.5) # TODO make radius dependent ?
        self.config["_buffer_size"] = buffer_size

        time_idx_start = 0
        time_idx_end = int((time_end-time_start)/time_step) + 1 + buffer_size
        time_ids = np.arange(time_idx_start, time_idx_end, 1)
        self.data["time"] = time_ids * self.time_step

    def _check_time_settings(self):
        # attributes
        time_config = self.config["model"]["time"]
        time_start = float(time_config["start"])
        time_end = float(time_config["end"])
        time_step = float(time_config["step"])
        time_block = float(time_config["block"])

        # check time settings
        if time_end <= time_start:
            raise ValueError("time_end should be larger than time_start.")
        if time_step < 0.0:
            raise ValueError("The timestep size can not be less than zero.")
        if modulo((time_end-time_start), time_block) > 1e-12:
            raise ValueError("Insertion interval (time_end - time_start) should be a multiple of t_unit.")
        if modulo(time_block, time_step) > 1e-12:
            raise ValueError("Variable t_unit should be a multiple of time_step.")

    def _check_mass_settings(self):
        # attributes
        mass_config = self.config["model"]["mass_g"]
        mass_target = float(mass_config["per_block"])
        mass_min = float(mass_config["min"])
        mass_max = float(mass_config["max"])
        mass_tol = float(mass_config["tol"])

        # check time settings
        if mass_target <= mass_min:
            raise ValueError("Target mass should be larger than minimum bubble mass")
        if mass_max <= mass_min:
            raise ValueError("Maximum bubble mass can not be larger than minimum bubble mass.")
        if mass_tol < 0.0:
            raise ValueError("Mass tolerance can not be less than zero.")

    def _check_modules(self):
        cfd_prog = str(self.config["packages"]["cfd_program"])
        if (cfd_prog.lower() != "openfoam") & (cfd_prog.lower() != "ansys_cfd"):
            raise RuntimeError("CFD program type not found. It should be either OpenFOAM or ANSYS_CFD")

    def _get_openfoam_type(self) -> str:
        """
        Obtain the OpenFOAM type by checking the cfd_version string according to predefined regex patterns.
        E.g.    cfd_version = v2312-foss-2023a => "com"
                cfd_version = 11-foss-2023a => "org"

        Args:
            cfd_version: The OpenFOAM version.

        Returns:
            Foundation ("com") or ESI ("org")

        Raises:
            ValueError: Raises an exception if the cfd_version is unknown.
        """

        cfd_version = self.config["packages"]["cfd_version"]

        com_pattern = re.compile("^(v[0-9]+-foss-[0-9]+[a-b])$")
        com_match = re.match(com_pattern, cfd_version)
        if com_match:
            return "com"

        org_pattern = re.compile("^([0-9]+).*(-foss-).*([0-9]+[a-b])$")
        org_match = re.match(org_pattern, cfd_version)
        if org_match:
            return "org"

        raise ValueError("The cfd_version variable is unknown for OpenFOAM. \
                            Use the format [v2312-foss-2023a] for com or [11-foss-2023a] for org version")

    def _initialize_cells(self):
        '''initialize arrays for time, alpha, and velocity'''
        n_faces = len(self.data["inlet_faces"])
        n_timesteps = len(self.data["time"])
        self.data["alpha"] = np.ones([n_faces, n_timesteps, 1], dtype=np.float64)
        self.data["velocity"] = np.ones([n_faces, n_timesteps, 3], dtype=np.float64)
        self.data["velocity"][:, :, :] *= self.velocity_bc * self.data["inlet_normal"][None, None, :]

    def _prepare_boundary_data_dir(self):
        cwd = self.cwd
        inlet_name = self.inlet_name
        inlet_faces = self.data["inlet_faces"]

        # Prepare OpenFOAM-directory
        self.logger.info("Creating 'boundaryData' folder.")
        constant_path = os.path.join(cwd, "constant")
        if not os.path.exists(constant_path):
            raise RuntimeError("'constant' folder does not exist in the OpenFOAM case directory.")
        self.config["_boundary_data_path"] = os.path.join(cwd, "constant", "boundaryData")

        # Create boundaryData
        os.mkdir(os.path.join(cwd, "constant", "boundaryData"))
        os.mkdir(os.path.join(cwd, "constant", "boundaryData", inlet_name))
        boundary_inlet_path = os.path.join(cwd, "constant", "boundaryData", inlet_name)
        # command = f"touch {boundary_inlet_path}/points"
        # subprocess.run(command, shell=True)

        # Write 'points'-file
        points_out_path = os.path.join(boundary_inlet_path, "points")
        n_faces = len(inlet_faces[:, 0])
        with open(points_out_path, 'w') as f:
            f.write(str(n_faces)+'\n')
            f.write('('+'\n')
            for fid in range(n_faces):
                f.write('(' + str(inlet_faces[fid, 1]) + ' ' + str(inlet_faces[fid, 2]) + ' ' + str(inlet_faces[fid, 3]) + ') \n')
            f.write(')'+'\n')


    def _iterate(self):
        '''model and write for every insertion block'''
        time_start = self.time_start
        time_end = self.time_end
        time_block = self.time_block
        n_blocks = int((time_end-time_start)/time_block)

        # initialize model buffers
        self.model.initialize_buffer()
        n_check = np.sum(self.model.buffer[:,:,0])
        temp_slice_t10 = self.data["alpha"][:,10,0]
        temp_slice_t11 = self.data["alpha"][:,11,0]

        self.logger.info(f"Discretized {n_blocks} intervals of {time_block} s between time_start {time_start} s and time_end {time_end} s.")
        for block_idx in range(n_blocks):
            n_alpha_old = np.sum(self.data["alpha"])
            # inlet modelling
            with self.logger.function_call(name="model"):
                self.model.initialize(block_idx)
                self.model.run(block_idx)

            # writer
            with self.logger.function_call(name="writer"):
                self.writer.initialize(block_idx)
                self.writer.run()

    def _save_csv(self):
        # print csv to visualize pre-inlet domain
        self.logger.info("Saving inlet profile to csv-files.")
        csv_file_path = os.path.join(self.config["_output_path"], "inlet_data.csv")
        inlet_all_variable = np.concatenate((self.data["alpha"][:,:,:], self.data["velocity"][:,:,:]), axis=2)
        csv_header = "x_coord,y_coord,z_coord,alpha,velocity_x,velocity_y,velocity_z"

        # get cell coordinates in x,y,z space
        n_faces = len(self.data["inlet_faces"])
        n_timesteps = len(self.data["time"])
        face_list_extended = np.array([self.data["inlet_faces"][:, 1:4]] * n_timesteps)
        time_velocity_product = np.tensordot(self.data["time"][:], self.velocity_bc * self.data["inlet_normal"][:], axes=0)
        cell_coords = face_list_extended[:, :, :] - time_velocity_product[:, None, : ]
        cell_coords = np.swapaxes(cell_coords, 0, 1)
        cell_coords = np.reshape(cell_coords, (n_timesteps * n_faces, -1), order='C')

        # save csv
        inlet_var_reshaped = np.reshape(inlet_all_variable, (n_timesteps * n_faces, -1), order="C")
        inlet_ds = np.concatenate((cell_coords, inlet_var_reshaped), axis=1)
        np.savetxt(csv_file_path, inlet_ds, fmt='%.6e',
                    header=csv_header, delimiter=",", comments='')
        self.logger.info("Inlet profile saved to csv-files.")
