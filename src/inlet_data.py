from utils import np
from utils import dataclass

@dataclass
class InletData:
    """Data storage class for SBM."""
    faces: np.array = None
    normal: np.array = None
    velocity: np.array = None
    alpha: np.array = None

    def __init__(self):
        self.clean_samples()

    def store_geometry(self, faces, normal):
        self.faces = faces
        self.normal = normal

    def initialize(self, config):
        '''initialize inlet data storage'''
        velocity = float(config["model"]["velocity"])
        timesteps_per_block = config["_timesteps_per_block"]
        buffer_size = config["_buffer_size"]
        n_time = timesteps_per_block + buffer_size

        n_faces = len(self.faces)
        self.alpha = np.ones([n_faces, n_time, 1], dtype=np.float64)
        self.velocity = np.ones([n_faces, n_time, 3], dtype=np.float64)
        self.velocity[:, :, :] *= velocity * self.normal[None, None, :]

    def store_block_data(self, alpha, velocity):
        self.alpha[:, :, :] = alpha[:, :, :]
        self.velocity[:, :, :] = velocity[:, :, :]

    def store_bubble(self, mass, position, is_defined=False):
        if is_defined:
            id = self._samples_moving_pointer
            samples = self.samples
        else:
            id = self._rejected_samples_moving_pointer
            samples = self.rejected_samples

        samples["mass"][id[0]] = mass
        samples["position"][id[0]] = position
        id[0] += 1

        # increase array sizes if half full
        if id[0] > len(samples["mass"])//2:
            samples["mass"] += [None]*len(samples["mass"])
            samples["position"] += [None]*len(samples["position"])

    def get_samples(self):
        '''Get samples and rejected samples dictionary with entries in numpy array format'''
        samples_end = self._samples_moving_pointer[0]
        rejected_samples_end = self._rejected_samples_moving_pointer[0]
        samples_clean, rejected_samples_clean = {}, {}
        for key in self.samples.keys():
            samples_clean[key] = self.samples[key][:samples_end]

        for key in self.rejected_samples.keys():
            rejected_samples_clean[key] = self.rejected_samples[key][:rejected_samples_end]

        return samples_clean, rejected_samples_clean

    def calc_geometry_vars(self):
        '''Calculate inlet bounds and center'''
        x = self.faces[:,1:4]
        self.min = np.min(x[:,:],axis=0)
        self.max = np.max(x[:,:],axis=0)
        self.center = np.average(x[:,:],axis=0)

    def clean_samples(self):
        # initialization of (mutable) dictionaries
        self.samples = {"mass": [None]*10000, "position": [None]*10000}
        self.rejected_samples = {"mass": [None]*10000, "position": [None]*10000}

        # pointers to array of samples and rejected samples
        self._samples_moving_pointer = [0]
        self._rejected_samples_moving_pointer = [0]
