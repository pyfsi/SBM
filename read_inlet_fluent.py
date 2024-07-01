import ansys.fluent.core as pyfluent
import sys
from pathlib import Path
from contextlib import redirect_stdout
import numpy as np
import multiprocessing

# Read input from bash-script
if len(sys.argv) != 9:
    sys.exit("8 arguments should be used: dimensions - case path - case name - start time step - number of time steps - "
             "inlet boundary name - max number of nodes per face - cores.")
dimensions = int(sys.argv[1])
case_path = Path(sys.argv[2])
case_name = str(sys.argv[3])
time_step_start = int(sys.argv[4])
n_time_steps = int(sys.argv[5])
boundary_name = str(sys.argv[6])
mnpf = int(sys.argv[7])
cores = min(int(sys.argv[8]), multiprocessing.cpu_count())

with open('read_inlet.log', 'w') as logfile:
    with redirect_stdout(logfile):
        # Launch fluent and read case
        session = pyfluent.launch_fluent(mode="solver", precision="double", version=f"{dimensions}d",
                                         processor_count=cores, show_gui=True)

        session.scheme_eval.scheme_eval("(enable-dynamic-mesh-node-ids #t)")
        session.tui.file.read_case(str(case_path / case_name))

        # Getting the thread id from Fluent similarly as is done in CoCoNuT
        # Could be replaced by using PyFluent SVARS, but does not work yet

        report = Path('report.sum')
        report.unlink(missing_ok=True)
        # session.tui.report.summary("yes", path.name)
        session.results.report.summary(write_to_file=True, file_name=report.name)
        boundary_names = {boundary_name: None}
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
                        raise RuntimeError(f'Format of {fp} not recognized')
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
        with open('udf_sbm_template.c', 'r') as infile:
            with open('udf_sbm.c', 'w') as outfile:
                for line in infile:
                    line = line.replace('|MAX_NODES_PER_FACE|', str(mnpf))
                    # one more than n_time_steps, because initial value is included
                    line = line.replace('|N_TIME_STEPS|', str(n_time_steps + 1))
                    line = line.replace('|TIME_STEP_START|', str(time_step_start))
                    outfile.write(line)

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
tol = 1e-14
if np.all(normal_inlets.std(axis=0) < tol):
    normal_inlet = normal_inlets.mean(axis=0)
else:
    raise NotImplementedError('SBM cannot deal with non-planar boundaries')

# Save inlet and normal in Python Numpy-array format
np.save("inletPython.npy", coord_list)
np.save("normalInletPython.npy", normal_inlet)
np.save("face_ids.npy", face_ids)
