import numpy as np
import matplotlib.pyplot as plt

R = 0.25  # m (radius of bubble)
D = 2*R  # m (diameter of bubble)
# Length of pre-domain in meter (normally calculated the other way around)
L = 4*R  # m
U = 0.5  # m/s (inlet velocity)
# Length of pre-domain in seconds (normally calculated the other way around)
T = L/U  # s 
dt = 1e-4  # s (time step of simulation)

n_slices = int(T/dt)  # number of mesh slices in pre-domain

# Read inlet discretization
faces = np.load('inletPython.npy')
normal = np.load('normalInletPython.npy')

n_faces = faces.shape[0]

VOFw = np.ones((n_faces, n_slices, 1))

midpoint = faces[:, 1:4].mean(axis=0)
midpoint += normal*U*T/2

for i in range(n_faces):
    for j in range(n_slices):
        point = faces[i, 1:4] + normal*U*dt*j
        distance = np.linalg.norm(midpoint-point)
        if distance < R:
            VOFw[i, j, 0] = 0


np.save('inletDefinition-VOFw.npy', VOFw)

# plot slice  # todo: tomography plotter
plt.figure()
plt.scatter(faces[:, 1], faces[:, 2], VOFw[:, 10000, 0])
