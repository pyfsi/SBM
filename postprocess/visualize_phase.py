from utils import plt, np, os
from utils import approx_equal

# !!! TODO: too heavy to visualize in Python.
# Import clean data and visualize in paraview
# in paraview: load .csv -> tables to point data -> !!!

def visualize_phase(coords, alpha, phase=0):
    x,y,z = coords
    is_phase = approx_equal(alpha, phase, tolerance=1e-3)

    temp = np.ones((3,3,3),
                    dtype=object) == 1

    colors = np.empty(temp.shape, dtype=object)
    colors[temp] = "blue"

    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(temp, facecolors=colors, edgecolor="k")
    ax.set(xlabel='Axis 1', ylabel='Axis 2', zlabel='Time')
    ax.set_aspect("auto")
    plt.show()

if __name__ == "__main__":
    # paths
    cwd = os.getcwd()
    time_path = os.path.join(cwd, "sbm_files", "inletDefinition-time.npy")
    inlet_path = os.path.join(cwd, "sbm_files", "inletPython.npy")
    alpha_path = os.path.join(cwd, "sbm_files", "inletDefinition-VOFw.npy")

    # read sbm vof
    time = np.load(time_path)
    inlet = np.load(inlet_path) # id - x - y - z - area
    inlet_coords = inlet[:,1:4] # x-y-z
    vof = np.load(alpha_path)[:, :-1] # cut last datapoint

    # PCA to get axes and non dimensional coords
    inlet_coords = inlet_coords.round(decimals=4)
    mean, var = np.mean(inlet_coords), np.var(inlet_coords)
    standard_coords = (inlet_coords - mean) / var
    covariance = np.cov(standard_coords.T)
    _, eigenvectors = np.linalg.eig(covariance) # NOTE: bad approximation

    # TODO
    new_coords = np.matmul(np.identity(3), standard_coords.T)
    x = np.unique(new_coords[0,:] * var + mean)
    y = np.unique(new_coords[1,:] * var + mean)
    grid = np.meshgrid(x, y, time)

    # remap VoF data according to grid shape
    vof_remapped = np.zeros(grid[0].shape)
    for i, xi in enumerate(x):
        for j, yi in enumerate(y):
            point_i = np.array(xi, yi)
            temp = inlet_coords[:,:-1] - point_i
            cell_id = np.argmin(np.linalg.norm(temp, ord=2, axis=1))
            vof_remapped[i,j,:] = vof[cell_id,:,0]

    # visualization
    visualize_phase(grid, vof_remapped)

    print("visualize_phase script ended successfully.")