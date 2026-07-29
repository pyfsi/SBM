import ansys.fluent.core as pyfluent
from utils import os, np, redirect_stdout, multiprocessing

def write_bc_ansyscfd(config):
    case_path = os.getcwd()
    dimensions = config["cfd"]["dimension"]
    case_name = str(config["fluent"]["case_name"])
    time_step_start = str(config["fluent"]["time_step_start"])
    cores = min(int(config["fluent"]["cores"]), multiprocessing.cpu_count())

    VOFw = np.load('inletDefinition-VOFw.npy')
    face_ids = np.load('face_ids.npy')
    np.savetxt(f'inlet_VOFw_start{time_step_start}.dat', VOFw.squeeze(), fmt='%i')
    np.savetxt('face_ids.dat', face_ids, fmt='%i', header=f'{face_ids.shape[0]}')

    log_path = os.path.join(os.getcwd(), "sbm_files", 'write_bc.log')
    with open(log_path, 'w') as logfile:
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
