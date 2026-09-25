from utils import np, hashlib

class SampleGenerator():
    '''
    Class to generate dispersed phase inclusions by modifying the values of a 3D grid.
    '''
    def __init__(self, config,):
        self.seed = str(config["model"].get("seed", None))
        self.rng = None
        self._initialize_rng()

    def sample_random(self, fid_bounds: list, tid_bounds: list, mass_bounds:list):
        # TODO change to output position and time instead of fid and tid
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

    # == Protected functions ==
    def _convert_string_to_int(self, a):
        '''String to int conversion method for seed generation based on Python's random package'''
        if isinstance(a, str):
            out = a.encode()
        else:
            out = a
        return int.from_bytes(out + hashlib.sha512(out).digest(), 'big')

    def _initialize_rng(self):
        # set seed for random number generator
        seed = self.seed
        if seed.lower()=="none":
            seed = None
        seed_msg = f"Using seed {seed} for random number generator"
        print(seed_msg)

        # set random number generator
        if isinstance(seed, str):
            seed_converted = self._convert_string_to_int(seed)
        self.rng = np.random.default_rng(seed=seed_converted)
