import numpy as np

VOFw = np.load('inletDefinition-VOFw.npy')
face_ids = np.load('face_ids.npy')
np.savetxt('inlet_VOFw.dat', VOFw.squeeze(), fmt='%i')
np.savetxt('face_ids.dat', face_ids, fmt='%i', header=f'{face_ids.shape[0]}')
