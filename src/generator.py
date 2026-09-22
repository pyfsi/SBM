from utils import np, os, random
from utils import PI

class Generator():
    '''
    Class to generate dispersed phase inclusions by modifying the values of a 3D grid.
    '''
    def __init__(self, config,):
        self.seed = str(config["model"].get("seed", None))
        self._set_seed()

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
        sample["face"] = random.randint(fid_bounds[0], fid_bounds[1])

        # sample time index for temporal coordinates
        sample["time"] = random.randint(tid_bounds[0], tid_bounds[1])

        # sample bubble mass
        sample["mass"] = random.uniform(mass_bounds[0], mass_bounds[1])

        return sample

    # == Protected functions ==
    def _set_seed(self):
        # set seed for random number generator
        seed = self.seed
        if seed.lower()=="none":
            seed = None
        seed_msg = f"Using seed {seed} for random number generator"
        print(seed_msg)
        random.seed(seed)
