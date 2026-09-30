from utils import os, shutil, re, np
from utils import modulo
from utils import PI

# SBM components
from .logger import Logger
from .preprocessor import Preprocessor
from .reader import Reader
from .model import Model
from .writer import Writer
from .plotter import Plotter
from .inlet_data import InletData

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
        self.density_gas = float(config["cfd"]["rho_g"])
        self.mg_max = float(config["model"]["mass_g"]["max"])

        self.purge_flag = bool(self.config["settings"]["purge_boundary_data"])
        self.plotter_flag = bool(self.config["settings"]["activate_plotter"])

        self.cwd = os.getcwd()

    def check_case(self):
        # check SBM config
        self._check_time_settings()
        self._check_mass_settings()
        self._check_modules()

    def purge_boundary_data(self):
        self.config["_boundary_data_path"] = os.path.join(self.cwd, "constant", "boundaryData")

        if self.purge_flag:
            if os.path.exists(self.config["_boundary_data_path"]):
                shutil.rmtree(self.config["_boundary_data_path"])
        else:
            raise RuntimeError("'boundaryData' directory already exists.")

    def initialize(self):
        self._make_output_dir()
        self._add_time_config()

        # === initialize sbm components ===
        # diagnostics
        self.logger = Logger(self.config)

        # inlet data storage
        self.inlet_data = InletData()

        # preprocessor
        self.preprocessor = Preprocessor(self.config)

        # reader
        self.reader = Reader(self.config, self.inlet_data, self.logger)

        # inlet modelling
        self.model = Model(self.config, self.inlet_data, self.logger)

        # writer
        self.writer = Writer(self.config, self.inlet_data, self.logger)
        self.writer.check()

        # plotter
        self.plotter = Plotter()

    def run(self):
        with self.logger.function_call(name="preprocessor"):
            self.preprocessor.initialize()
            self.preprocessor.run()

        with self.logger.function_call(name="reader"):
            self.reader.initialize()
            self.reader.run()

        # initialize inlet data storage
        self.inlet_data.initialize(self.config)
        self.inlet_data.calc_geometry_vars()

        # prepare directory for boundary condition
        self._prepare_boundary_data_dir()

        # main SBM iteration
        self._iterate()

        # save plot
        if self.plotter_flag:
            self.plotter.save_plot()

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
        cwd = self.cwd
        sbm_output_path = os.path.join(cwd, "sbm_files")
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

    def _prepare_boundary_data_dir(self):
        cwd = self.cwd
        inlet_name = self.inlet_name
        inlet_faces = self.inlet_data.faces

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

        # singular intialization routines
        self.model.initialize_buffer()
        inlet_data = self.inlet_data
        self.plotter.initialize(self.config, inlet_data)

        self.logger.info(f"Discretized {n_blocks} intervals of {time_block} s between time_start {time_start} s and time_end {time_end} s.")
        for block_idx in range(n_blocks):
            # inlet modelling
            with self.logger.function_call(name="model"):
                self.model.initialize(block_idx)
                self.model.run()

            # writer
            with self.logger.function_call(name="writer"):
                self.writer.initialize(block_idx)
                self.writer.run()

            # plotter
            if self.plotter_flag:
                with self.logger.function_call(name="plotter"):
                    samples, rejected_samples = self.inlet_data.get_samples()
                    self.plotter.plot(samples, rejected_samples)

            # inlet data cleanup
            self.inlet_data.clean_samples()

    def _add_time_config(self):
        time_step = self.time_step
        time_block = self.time_block
        mg_max = self.mg_max
        density_gas = self.density_gas
        velocity_bc = self.velocity_bc

        # calculate bufer size based on radius
        radius_bubble = ((3.0*mg_max)/(4.0*PI*density_gas))**(1.0/3.0)
        buffer_size = int((radius_bubble/velocity_bc)/time_step)

        timesteps_per_block = int(time_block/time_step)
        self.config["_timesteps_per_block"] = timesteps_per_block
        buffer_size = int(timesteps_per_block + buffer_size)
        self.config["_buffer_size"] = buffer_size

    # === check functions ===
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
