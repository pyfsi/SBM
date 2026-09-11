from utils import np
from utils import dataclass

@dataclass
class InletData:
    """Data storage class for SBM."""
    faces: np.array = None
    normal: np.array = None
    velocity: np.array = None
    alpha: np.array = None

    def store_geometry(self, faces, normal):
        self.faces = faces
        self.normal = normal

    def store_block_data(self, alpha, velocity):
        self.alpha[:, :, :] = alpha[:, :, :]
        self.velocity[:, :, :] = velocity[:, :, :]

    def initialize_block_data(self, config):
        '''initialize arrays for alpha and velocity'''
        velocity = float(config["model"]["velocity"])
        timesteps_per_block = config["_timesteps_per_block"]
        buffer_size = config["_buffer_size"]
        n_time = timesteps_per_block + buffer_size

        n_faces = len(self.faces)
        self.alpha = np.ones([n_faces, n_time, 1], dtype=np.float64)
        self.velocity = np.ones([n_faces, n_time, 3], dtype=np.float64)
        self.velocity[:, :, :] *= velocity * self.normal[None, None, :]