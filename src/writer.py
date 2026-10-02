from utils import os, np, pd, subprocess, partial, ThreadPoolExecutor

class Writer():
    def __init__(self, config, data, logger):
        # configuration
        self.density_gas = float(config["cfd"]["rho_g"])
        self.time_start = float(config["model"]["time"]["start"])
        self.time_step = float(config["model"]["time"]["step"])
        self.inlet_name = str(config["cfd"]["inlet_name"])
        self.alpha_name = "alpha."+config["cfd"]["alpha_name"]
        self.velocity_bc = float(config["model"]["velocity"])
        self.save_csv = bool(config["settings"].get("save_csv", False))

        # paths
        self.cwd = os.getcwd()
        self.output_path = str(config.get("_output_path"))
        self.boundary_data_path = str(config.get("_boundary_data_path"))

        self.timesteps_per_block = int(config.get("_timesteps_per_block"))
        self.buffer_size = int(config.get("_buffer_size"))

        # reference to data storage
        self.inlet_data = data

        # logger
        self.logger = logger

    def initialize(self, block_idx):
        # absolute time index and time array
        abs_time_idx_block_start = self.timesteps_per_block * block_idx
        abs_time_idx_block_end = self.timesteps_per_block * (block_idx+1)
        abs_time_idx = np.arange(abs_time_idx_block_start, abs_time_idx_block_end, 1)
        self.time = abs_time_idx * self.time_step + self.time_start
        self.block_idx = block_idx

        self.inlet_faces = self.inlet_data.faces[:,:]
        self.inlet_normal = self.inlet_data.normal[:]
        self.alpha = self.inlet_data.alpha[:,:,:]
        self.velocity = self.inlet_data.velocity[:,:,:]

    def run(self):
        self._write_boundary_data()

        if self.save_csv:
            self._save_csv()

    def run_last(self, block_idx):
        ''''Write boundary condition for last time step'''
        # absolute time index and time array
        abs_time_idx_block_start = self.timesteps_per_block * block_idx
        abs_time_idx_block_end = self.timesteps_per_block * (block_idx+1)
        abs_time_idx = np.arange(abs_time_idx_block_start, abs_time_idx_block_end, 1)
        self.time = abs_time_idx * self.time_step + self.time_start
        self.block_idx = block_idx

        # set inlet data
        self.inlet_faces = self.inlet_data.faces[:,:]
        self.inlet_normal = self.inlet_data.normal[:]
        self.alpha = np.ones(self.inlet_data.alpha[:,:,:].shape)
        vel_bc = self.inlet_normal[:] * self.velocity_bc
        self.velocity = np.ones(self.inlet_data.velocity.shape[:2])[:,:,None] * vel_bc

        self._write_boundary_data()

    def check(self):
        '''Check boundary condition definition.'''
        cwd = os.getcwd()
        time_start = f"{self.time_start:g}"
        inlet_name = self.inlet_name
        alpha_name = self.alpha_name

        # Check boundary condition for 'U'. Boundary condition should be type timeVaryingMappedFixedValue
        command = f"grep -nrw {inlet_name} {time_start}/U | cut -d : -f 1 > lineNr_U"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNrU_path = os.path.join(cwd, "lineNr_U")
        lineNameNr = int(open(lineNrU_path, 'r').readline())
        lineTypeU = lineNameNr+1  # inlet BC type is defined on this line
        U_path = os.path.join(cwd, time_start, "U")
        readTypeU = (open(U_path).readlines())[lineTypeU]
        boundaryConditionU = (readTypeU.split())[-1][0:-1]
        os.remove(lineNrU_path)
        if boundaryConditionU != "timeVaryingMappedFixedValue":
            raise RuntimeError(f"The condition for 'U' at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \n")

        # Boundary condition should not be setAverage
        command = f"grep -nr 'setAverage' {time_start}/U | cut -d : -f 1 > lineNr_setAvg"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNrsetAvg_path = os.path.join(cwd, "lineNr_setAvg")
        try:
            lineNameNr = int(open(lineNrsetAvg_path, 'r').readline())-1  # Line where 'setAverage' is defined
            readSetAvg = (open(U_path).readlines())[lineNameNr]
            setAvg = (readSetAvg.split())[-1][0:-1]
            os.remove(lineNrsetAvg_path)
            if setAvg != "false":
                raise RuntimeError(f"The boundary condition at '{inlet_name}' defines an averaging operation for 'U'. \
                        This is not compatible with the transient inlet modelling defined in the Python script. \n")
        except:
            os.remove(lineNrsetAvg_path)

        # Check boundary condition for alpha
        command = f"grep -nrw {inlet_name} {time_start}/{alpha_name} | cut -d : -f 1 > lineNr_VOFw"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNr_VOFw_path = os.path.join(cwd, "lineNr_VOFw")
        lineNameNr = int(open(lineNr_VOFw_path, 'r').readline())
        lineTypeVOFw = lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
        alpha_path = os.path.join(cwd, str(time_start), alpha_name)
        readTypeVOFw = (open(alpha_path).readlines())[lineTypeVOFw]
        boundaryConditionVOFw = (readTypeVOFw.split())[-1][0:-1]
        os.remove(lineNr_VOFw_path)
        if boundaryConditionVOFw != "timeVaryingMappedFixedValue":
            raise RuntimeError(f"The condition for {alpha_name} at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \
                    See line {lineNameNr} \n")
        command = f"grep -nr 'setAverage' {time_start}/{alpha_name} | cut -d : -f 1 > lineNr_setAvg"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNr_setAvg_path = os.path.join(cwd, "lineNr_setAvg")
        try:
            lineNameNr = int(lineNr_setAvg_path)-1 #Line where 'setAverage' is defined
            readSetAvg = (open(alpha_path).readlines())[lineNameNr]
            setAvg = (readSetAvg.split())[-1][0:-1]
            os.remove(lineNr_setAvg_path)
            if setAvg != "false":
                raise RuntimeError(f"Error! The boundary condition at '{inlet_name}' defines an averaging operation for {alpha_name}. \
                        This is not compatible with the transient inlet modelling defined in the Python script.")
            print(f"Inlet definition in folder {time_start} is OK.")
        except:
            os.remove(lineNr_setAvg_path)

    # === Protected functions ===
    def _write_boundary_data(self):
        # attributes
        boundary_data_path = self.boundary_data_path
        boundary_inlet_path = os.path.join(boundary_data_path, "inlet")
        alpha_name = self.alpha_name
        n_faces = len(self.inlet_faces[:, 0])
        velocity = self.velocity
        alpha = self.alpha
        time = self.time

        # Write velocity and volume fraction for water at each time step.
        self.logger.info("Writing boundary condition to folder 'boundaryData'.")
        args = [boundary_inlet_path, alpha_name, n_faces, velocity, alpha, time]

        # with partial
        partial_write = partial(self.__write_time_step, args=args)
        n_workers = max(8, os.cpu_count())
        self.logger.info(f"Using {n_workers} threads for writing routine.")
        time_steps = [tid for tid in range(len(time))]
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            executor.map(partial_write, time_steps)

        self.logger.info("Boundary condition was successfully saved in 'boundaryData'.")

    def _save_csv(self):
        # number of face and timesteps
        n_faces = len(self.inlet_faces)
        n_timesteps = len(self.time)

        # print csv to visualize pre-inlet domain
        self.logger.info("Saving inlet profile to csv-files.")
        csv_name = f"inlet_data_block{self.block_idx}.csv"
        csv_file_path = os.path.join(self.output_path, csv_name)
        inlet_all_variable = np.concatenate((self.alpha[:,:n_timesteps,:], self.velocity[:,:n_timesteps,:]), axis=2)
        csv_header = "x_coord,y_coord,z_coord,alpha,velocity_x,velocity_y,velocity_z"

        # get cell coordinates in x,y,z space
        face_list_extended = np.array([self.inlet_faces[:, 1:4]] * n_timesteps)
        time_velocity_product = np.tensordot(self.time[:], self.velocity_bc * self.inlet_normal[:], axes=0)
        cell_coords = face_list_extended[:, :, :] - time_velocity_product[:, None, : ]
        cell_coords = np.swapaxes(cell_coords, 0, 1)
        cell_coords = np.reshape(cell_coords, (n_timesteps * n_faces, -1), order='C')

        # save csv
        inlet_var_reshaped = np.reshape(inlet_all_variable, (n_timesteps * n_faces, -1), order="C")
        inlet_ds = np.concatenate((cell_coords, inlet_var_reshaped), axis=1)
        np.savetxt(csv_file_path, inlet_ds, fmt='%.6e',
                    header=csv_header, delimiter=",", comments='')
        self.logger.info("Inlet profile saved to csv-files.")

    # === Private functions ===
    def __write_time_step(self, tid, args):
        # arguments
        boundary_inlet_path = args[0]
        alpha_name = args[1]
        n_faces = args[2]
        velocity = args[3]
        alpha = args[4]
        time = args[5]
        time_i = time[tid]
        trunc_time_i = f"{time_i:g}"

        # Set location of files
        U_out_path = f"{boundary_inlet_path}/{trunc_time_i}/U"
        VOFw_out_path = f"{boundary_inlet_path}/{trunc_time_i}/{alpha_name}"

        # make directory for each time step
        time_i_dir = os.path.join(boundary_inlet_path, str(trunc_time_i))
        os.mkdir(time_i_dir)
        make_u_dir = f"touch {boundary_inlet_path}/{trunc_time_i}/U;"
        subprocess.run(make_u_dir, shell=True)
        make_alpha_dir = f"touch {boundary_inlet_path}/{trunc_time_i}/{alpha_name};"
        subprocess.run(make_alpha_dir, shell=True)

        # write content
        with open(U_out_path, 'a+') as f:
            f.write('//Data points'+'\n')
            f.write(str(n_faces)+'\n'+'('+'\n')
            for fid in range(n_faces):
                f.write('(' + str(velocity[fid, tid, 0]) + ' ' + str(velocity[fid, tid, 1]) + ' ' + str(velocity[fid, tid, 2]) + ') \n')
            f.write(')'+'\n')

        with open(VOFw_out_path, 'a+') as f:
            f.write('//Data points'+'\n')
            f.write(str(n_faces)+'\n'+'('+'\n')
            for fid in range(n_faces):
                f.write(str(alpha[fid,tid,0])+'\n')
            f.write(')'+'\n')
