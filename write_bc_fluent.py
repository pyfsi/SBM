import ansys.fluent.core as pyfluent
import sys
from pathlib import Path
from contextlib import redirect_stdout
import numpy as np
import multiprocessing

# Read input from bash-script
if len(sys.argv) != 7:
    sys.exit("6 arguments should be used: dimensions - case path - case name - start time step - inlet boundary name"
             " - cores.")
dimensions = int(sys.argv[1])
case_path = Path(sys.argv[2])
case_name = str(sys.argv[3])
time_step_start = int(sys.argv[4])
boundary_name = str(sys.argv[5])
cores = min(int(sys.argv[6]), multiprocessing.cpu_count())

VOFw = np.load('inletDefinition-VOFw.npy')
face_ids = np.load('face_ids.npy')
np.savetxt(f'inlet_VOFw_start{time_step_start}.dat', VOFw.squeeze(), fmt='%i')
np.savetxt('face_ids.dat', face_ids, fmt='%i', header=f'{face_ids.shape[0]}')

with open('write_bc.log', 'w') as logfile:
    with redirect_stdout(logfile):
        # Launch fluent and read case
        session = pyfluent.launch_fluent(mode="solver", precision="double", version=f"{dimensions}d",
                                         processor_count=cores, show_gui=True, cleanup_on_exit=False)

        session.scheme_eval.scheme_eval("(enable-dynamic-mesh-node-ids #t)")
        if time_step_start:
            session.file.read(file_type="case-data", file_name=case_path / case_name)
        else:
            session.file.read(file_type="case", file_name=case_path / case_name)

        # Compile and use UDF
        session.tui.define.user_defined.compiled_functions("compile", "sbm_lib", "yes", "udf_sbm.c")
        session.tui.define.user_defined.compiled_functions("load", "sbm_lib")
        session.tui.define.user_defined.execute_on_demand('"read_vof_face_ids::sbm_lib"')

        # session.tui.define.boundary_conditions.set.velocity_inlet("inlet", (), "water", "volume-fraction", "yes",
        #                                                           "yes", "udf", "sbm_profile::libudf", "quit")
        session.setup.boundary_conditions.velocity_inlet['inlet'].phase['water'].volume_fraction = \
            {"option": "udf", "udf": "sbm_profile::sbm_lib"}

        session.solution.run_calculation.transient_controls.time_step_count = VOFw.shape[1] - 1
        session.exit()
