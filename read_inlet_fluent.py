import ansys.fluent.core as pyfluent
import numpy as np
import sys
from pathlib import Path

# Read input from bash-script
if len(sys.argv) != 6:
    sys.exit("5 arguments should be used: dimensions - case path - start time step - number of time steps - max number"
             "of nodes per face.")
dimensions = int(sys.argv[1])
case_path = Path(sys.argv[2])
time_step_start = int(sys.argv[3])
n_time_steps = int(sys.argv[4])
mnpf = int(sys.argv[5])

# Launch fluent and read case
session = pyfluent.launch_fluent(mode="solver", precision="double", version=f"{dimensions}d", processor_count=4,
                                 show_gui=True)

session.scheme_eval.scheme_eval("(enable-dynamic-mesh-node-ids #t)")
session.read_case(case_path)

# Getting the thread id from Fluent similarly as is done in CoCoNuT
# Could be replaced by using PyFluent SVARS, but does not work yet

report = Path('report.sum')
report.unlink(missing_ok=True)
# session.tui.report.summary("yes", path.name)
session.results.report.summary(write_to_file=True, file_name=report.name)
boundary_names = {'inlet': None}
check = 0
names_found = []
with open(report, 'r') as fp:
    for line in fp:
        if check == 3 and line.islower():
            line_list = line.strip().split()
            if len(line_list) == 3:
                name, thread_id, _ = line_list
            elif len(line_list) == 4:
                name, _, thread_id, _ = line_list
            else:
                raise RuntimeError(f'Format of {file} not recognize')
            if name in boundary_names and name not in names_found:
                boundary_names['inlet'] = thread_id
                names_found.append(name)
        if check == 3 and not line.islower():
            break
        if check == 2:  # skip 1 line
            check = 3
        if 'name' in line and check == 1:
            check = 2
        if 'Boundary Conditions' in line:
            check = 1
with open('inlets.txt', 'w') as fp:
    fp.write(f'{len(names_found)}\n')
    for name, thread_id in boundary_names.items():
        fp.write(f'{name} {thread_id}\n')

# Make substitutions in UDF
with open('udf_sbm.c') as infile:
    with open(case_path / 'udf_sbm.c') as outfile:
        for line in infile:
            line = line.replace('|MAX_NODES_PER_FACE|', str(mnpf))
            line = line.replace('|N_TIME_STEPS|', str(n_time_steps))
            line = line.replace('|TIME_STEP_START|', str(time_step_start))

# Compile and use UDF
session.tui.define.user_defined.compiled_functions("compile", "sbm_lib", "yes", "udf_sbm.c")
session.tui.define.user_defined.compiled_functions("load", "sbm_lib")
session.tui.define.user_defined.execute_on_demand('"get_inlet_thread_ids::sbm_lib"')
session.tui.define.user_defined.execute_on_demand('"store_faces_normals_ids::sbm_lib"')

# Close Fluent
session.exit()

# Read in faces file
faces = np.loadtxt('faces.dat', skiprows=1)
faces_n = faces.shape[0]  # number of faces is the number of rows
coord_list = np.ones((faces_n, 5))  # initialize on 1
coord_list[:, 0] = np.arange(faces_n)
coord_list[:, 1:4] = faces[:, 0:3]
coord_list[:, 4] = np.linalg.norm(faces[:, 3:6], axis=1)  # face areas
face_ids = faces[:, 6:10]

normal_inlets = faces[:, 3:6] / coord_list[:, 4].reshape(faces_n, 1)  # compute unit normals on each face
tol = 1e-15
if np.all(normal_inlets.std(axis=0) < tol):
    normal_inlet = normal_inlets.mean(axis=0)
else:
    raise NotImplementedError('SBM cannot deal with non-planar boundaries')

# Save inlet and normal in Python Numpy-array format
np.save("inletPython.npy", coord_list)
np.save("normalInletPython.npy", normal_inlet)
np.save("face_ids.npy", face_ids)
