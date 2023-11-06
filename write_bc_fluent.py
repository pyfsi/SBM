import ansys.fluent.core as pyfluent
import sys
import numpy as np
from pathlib import Path

# Read input from bash-script
if len(sys.argv) != 4:
    sys.exit("3 arguments should be used: dimensions - case path - start time step.")
dimensions = int(sys.argv[1])
case_path = Path(sys.argv[2])
time_step_start = int(sys.argv[3])


VOFw = np.load('inletDefinition-VOFw.npy')
face_ids = np.load('face_ids.npy')
np.savetxt(f'inlet_VOFw_start{time_step_start}.dat', VOFw.squeeze(), fmt='%i')
np.savetxt('face_ids.dat', face_ids, fmt='%i', header=f'{face_ids.shape[0]}')

# Launch fluent and read case
session = pyfluent.launch_fluent(mode="solver", precision="double", version=f"{dimensions}d", processor_count=4,
                                 show_gui=True)

session.scheme_eval.scheme_eval("(enable-dynamic-mesh-node-ids #t)")
session.read_case(case_path)
session.tui.define.user_defined.execute_on_demand('"read_vof_face_ids::sbm_lib"')
