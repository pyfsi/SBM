from utils import os, shutil, re
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

        # make sbm output dir
        self.cwd = os.getcwd()
        sbm_output_path = os.path.join(self.cwd, "sbm_files")
        # remove previous sbm output files
        if os.path.exists(sbm_output_path):
            shutil.rmtree(sbm_output_path)
            # print(f"Deleting previous SBM output files.")
        os.mkdir(sbm_output_path)

        # === add attributes to config ===
        self.config["_output_path"] = sbm_output_path
        is_cfd_solver_openfoam = (self.config["packages"]["cfd_program"].lower() == "openfoam")
        if is_cfd_solver_openfoam:
            self.config["_openfoam_type"] = self._get_openfoam_type()

        # === data storage ===
        self.data = {}
        self.data["inlet_faces"] = None
        self.data["inlet_normal"] = None
        self.data["time"] = None
        self.data["alpha"] = None
        self.data["velocity"] = None

    def check_settings(self):
        self._check_time_settings()
        self._check_mass_settings()
        self._check_modules()

    def purge_previous(self):
        purge_boundary_data = self.config["settings"]["purge_boundary_data"]
        if purge_boundary_data:
            boundary_data_path = os.path.join(self.cwd, "constant", "boundaryData")
            if os.path.exists(boundary_data_path):
                shutil.rmtree(boundary_data_path)
        else:
            raise RuntimeError("'boundaryData' directory already exists.")

    def initialize(self):
        '''Initialize SBM components.'''
        config = self.config
        data = self.data

        # === initialize all sbm components ===
        # diagnostics
        # self.profiler = Profiler(config)
        self.logger = Logger(config)

        # preprocessor
        self.preprocessor = Preprocessor(config)

        # reader
        self.reader = Reader(config, data)

        # inlet modelling
        self.model = Model(config, data, self.logger)

        # writer
        self.writer = Writer(config, data, self.logger)

        # plotter
        self.plotter = None

    def run(self):
        with self.logger.function_call(name="preprocessor"):
            self.preprocessor.run()

        with self.logger.function_call(name="reader"):
            self.reader.initialize()
            self.reader.run()

        with self.logger.function_call(name="model"):
            self.model.initialize()
            self.model.run()

        with self.logger.function_call(name="writer"):
            self.writer.initialize()
            self.writer.run()

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

        for file in filenames:
            # move area postProcess file
            target_area_path = os.path.join(output_path, file)
            if os.path.exists(target_area_path):
                os.remove(target_area_path)
            area_path = os.path.join(cwd, "0", file)
            if os.path.exists(area_path):
                shutil.move(area_path, output_path)

    # == protected functions ==
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

