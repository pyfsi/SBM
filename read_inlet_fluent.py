import ansys.fluent.core as pyfluent

session = pyfluent.launch_fluent(mode="solver", precision="double", version="3d", processor_count=4, show_gui=True)

session.scheme_eval.scheme_eval("(enable-dynamic-mesh-node-ids #t)")
session.file.read_case_data(file_type="case-data", file_name="first_test")


# Getting the thread id from Fluent similarly as is done in CoCoNuT
# Could be replaced by using PyFluent SVARS, but does not work yet
file = 'report.sum'
session.results.report.summary(write_to_file=True, file_name="report.sum")
boundary_names = {'inlet': None}
check = 0
names_found = []
with open(file, 'r') as fp:
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
			check =2
		if 'Boundary Conditions' in line:
			check = 1
with open('inlets.txt', 'w') as fp:		 
	fp.write(f'{len(names_found)}\n')
	for name, id in boundary_names.items():
		fp.write(f'{name} {id}\n')

session.tui.define.user_defined.compiled_functions("compile", "libudf","yes", "udf_sbm.c")
session.tui.define.user_defined.compiled_functions("load", "libudf")
session.tui.define.user_defined.execute_on_demand('"get_inlet_thread_ids::libudf"')
session.tui.define.user_defined.execute_on_demand('"store_faces_normals_ids::libudf"')
