from utils import np, hashlib

class SampleGenerator():
    '''
    Class to generate dispersed phase inclusions by modifying the values of a 3D grid.
    '''
    def __init__(self, config,):
        self.seed = str(config["model"].get("seed", None))
        self.rng = self._initialize_rng()
        self.sampling_config = config["model"]["sampling"]
        self.sampling_pdfs = self._initialize_pdfs()

    def sample_random_discrete(self, fid_bounds: list, tid_bounds: list, mass_bounds:list):
        '''
        Sample bubble parameters from predefined statistical distributions.
        Args:
            fid_bounds: lower and upper bounds for face index
            tid_bounds: lower and upper bounds for time index
            mass_bounds: lower and upper bounds for mass of one bubble

        Returns:
            Returns face_idx, time_idx, mass
        '''
        sample = {"face": None, "time": None, "mass": None}

        # sample cell index for spatial coordinates
        sample["face"] = self.rng.integers(low=fid_bounds[0], high=fid_bounds[1])

        # sample time index for temporal coordinates
        sample["time"] = self.rng.integers(low=tid_bounds[0], high=tid_bounds[1])

        # sample bubble mass
        sample["mass"] = self.rng.uniform(low=mass_bounds[0], high=mass_bounds[1])

        return sample

    def sample_random(self, position_bounds: list, time_bounds: list, mass_bounds:list):
        '''
        Sample bubble parameters from predefined statistical distributions.
        Args:
            position_bounds: lower and upper bounds for position
            time_bounds: lower and upper bounds for time
            mass_bounds: lower and upper bounds for mass of one bubble

        Returns:
            Returns sample variables defined by position, time, and mass
        '''
        sample = {"position": None, "mass": None}

        # sample position
        pos_pdf = self.sampling_pdfs["position"]["pdf"]
        pos_pdf_params = self.sampling_pdfs["position"]["params"]
        pos_sample = pos_pdf(**pos_pdf_params, size=3)
        pos_sample_corr = pos_sample * (position_bounds[1]-position_bounds[0]) + position_bounds[0]
        sample["position"] = self.__clamp(pos_sample_corr, position_bounds[0], position_bounds[1])

        # sample time
        time_pdf = self.sampling_pdfs["time"]["pdf"]
        time_pdf_params = self.sampling_pdfs["time"]["params"]
        time_sample = time_pdf(**time_pdf_params)
        time_sample_corr = time_sample * (time_bounds[1]-time_bounds[0]) + time_bounds[0]
        sample["time"] = self.__clamp(time_sample_corr, time_bounds[0], time_bounds[1])

        # sample bubble mass
        mass_pdf = self.sampling_pdfs["mass"]["pdf"]
        mass_pdf_params = self.sampling_pdfs["mass"]["params"]
        mass_sample = mass_pdf(**mass_pdf_params)
        mass_sample_corr = mass_sample * (mass_bounds[1]-mass_bounds[0]) + mass_bounds[0]
        sample["mass"] = self.__clamp(mass_sample_corr, mass_bounds[0], mass_bounds[1])

        return sample

    # == Private functions ==
    @staticmethod
    def __clamp(val, lb, ub):
        return np.minimum(np.maximum(val, lb), ub)

    @staticmethod
    def __convert_string_to_int(a):
        '''String to int conversion method for seed generation based on Python's random package'''
        if isinstance(a, str):
            out = a.encode()
        else:
            out = a
        return int.from_bytes(out + hashlib.sha512(out).digest(), 'big')

    # == Protected functions ==
    def _initialize_rng(self):
        # set seed for random number generator
        seed = self.seed
        if seed.lower()=="none":
            seed = None
        seed_msg = f"Using seed {seed} for random number generator"
        print(seed_msg)

        # set random number generator
        if isinstance(seed, str):
            seed_converted = self.__convert_string_to_int(seed)

        return np.random.default_rng(seed=seed_converted)

    def _initialize_pdfs(self):
        # init pdf output
        sampling_pdfs = {
            "position": {"pdf":None, "params":{}},
            "time": {"pdf":None, "params":{}},
            "mass": {"pdf":None, "params":{}},
        }

        # define mapping between string to distribution function
        pdf_dict = {
            "uniform": self.rng.uniform,
            "normal": self.rng.normal,
            "beta": self.rng.beta,
        }

        # set distribution types and params
        for var in sampling_pdfs.keys():
            pdf_type = str(self.sampling_config[var].get("type", "uniform"))
            if pdf_type not in pdf_dict.keys():
                raise RuntimeWarning(f"Probability distribution of type {pdf_type} not defined. Using uniform distribution.")

            # type
            sampling_pdfs[var]["pdf"] = pdf_dict[pdf_type]

            # params
            match pdf_type:
                case "uniform":
                    # default := bounded by [0,1] with equal probability for all values
                    sampling_pdfs[var]["params"]["low"] = 0.0
                    sampling_pdfs[var]["params"]["high"] = 1.0
                case "normal":
                    # default := unbounded, values outside of bounds [0,1] very unlikely
                    sampling_pdfs[var]["params"]["loc"] = float(self.sampling_config[var].get("loc", 0.5))
                    sampling_pdfs[var]["params"]["scale"] = float(self.sampling_config[var].get("scale", 0.15))
                case "beta":
                    # default := bounded by [0,1], peak at 0.5
                    sampling_pdfs[var]["params"]["a"] = float(self.sampling_config[var].get("a", 2.0))
                    sampling_pdfs[var]["params"]["b"] = float(self.sampling_config[var].get("b", 2.0))

        return sampling_pdfs

