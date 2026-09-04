from utils import os, subprocess, shutil, partial, ThreadPoolExecutor
from utils import truncate

class Writer():
    def __init__(self, config, data, logger):
        # configuration
        self.density_gas = float(config["cfd"]["rho_g"])
        self.time_start = str(config["model"]["time"]["start"])
        self.time_step = float(config["model"]["time"]["step"])
        self.inlet_name = str(config["cfd"]["inlet_name"])
        self.alpha_name = "alpha."+config["cfd"]["alpha_name"]

        # paths
        self.cwd = os.getcwd()
        self.output_path = str(config.get("_output_path"))

        # reference to data storage
        self.data = data
        self.face_list = data["inlet_faces"]

        # logger
        self.logger = logger

    def initialize(self):
        pass

    def run(self):
        self._prepare_boundary_data_dir()
        self._write_boundary_data()

    def check(self):
        '''Check boundary condition definition.'''
        cwd = os.getcwd()
        time_start = self.time_start
        inlet_name = self.inlet_name
        alpha_name = self.alpha_name

        # Check boundary condition for 'U'. Boundary condition should be type timeVaryingMappedFixedValue
        command = f"grep -nrw {inlet_name} {time_start}/U | cut -d : -f 1 > lineNr_U"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNrU_path = os.path.join(cwd, "lineNr_U")
        lineNameNr = int(open(lineNrU_path, 'r').readline())
        lineTypeU = lineNameNr+1  # inlet BC type is defined on this line
        U_path = os.path.join(cwd, str(time_start), "U")
        readTypeU = (open(U_path).readlines())[lineTypeU]
        boundaryConditionU = (readTypeU.split())[-1][0:-1]
        os.remove(lineNrU_path)
        if boundaryConditionU != "timeVaryingMappedFixedValue":
            raise RuntimeError(f"The condition for 'U' at the boundary '{inlet_name}' is not set to 'timeVaryingMappedFixedValue'. \n")

        # Boundary condition should not be setAverage
        command = f"grep -nr 'setAverage' {time_start}/U | cut -d : -f 1 > lineNr_setAvg"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNrsetAvg_path = os.path.join(cwd, "lineNr_setAvg")
        lineNameNr = int(open(lineNrsetAvg_path, 'r').readline())-1  # Line where 'setAverage' is defined
        readSetAvg = (open(U_path).readlines())[lineNameNr]
        setAvg = (readSetAvg.split())[-1][0:-1]
        os.remove(lineNrsetAvg_path)
        if setAvg != "false":
            raise RuntimeError(f"The boundary condition at '{inlet_name}' defines an averaging operation for 'U'. \
                    This is not compatible with the transient inlet modelling defined in the Python script. \n")

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
    def _prepare_boundary_data_dir(self):
        cwd = self.cwd
        inlet_name = self.inlet_name
        face_list = self.data["inlet_faces"]

        # Prepare OpenFOAM-directory
        self.logger.info("Creating 'boundaryData' folder.")
        constant_path = os.path.join(cwd, "constant")
        if not os.path.exists(constant_path):
            raise RuntimeError("'constant' folder does not exist in the OpenFOAM case directory.")
        self.boundary_data_path = os.path.join(cwd, "constant", "boundaryData")
        if os.path.exists(self.boundary_data_path):
            raise RuntimeError("'boundaryData' folder already exists.")

        # Create boundaryData
        os.mkdir(os.path.join(cwd, "constant", "boundaryData"))
        os.mkdir(os.path.join(cwd, "constant", "boundaryData", inlet_name))
        boundary_inlet_path = os.path.join(cwd, "constant", "boundaryData", inlet_name)
        command = f"touch {boundary_inlet_path}/points"
        subprocess.run(command, shell=True)

        # Write 'points'-file
        points_out_path = os.path.join(boundary_inlet_path, "points")
        n_faces = len(face_list[:, 0])
        with open(points_out_path, 'a+') as f:
            f.write(str(n_faces)+'\n')
            f.write('('+'\n')
            for j in range(n_faces):
                f.write('(' + str(face_list[j, 1]) + ' ' + str(face_list[j, 2]) + ' ' + str(face_list[j, 3]) + ') \n')
            f.write(')'+'\n')

    def _write_boundary_data(self):
        # attributes
        time_step = self.time_step
        boundary_inlet_path = os.path.join(self.boundary_data_path, "inlet")
        alpha_name = self.alpha_name
        n_faces = len(self.data["inlet_faces"][:, 0])
        U_inlet = self.data["velocity"]
        VOF_inlet = self.data["alpha"]
        time_inlet = self.data["time"]

        # Write velocity and volume fraction for water at each time step.
        self.logger.info("Writing boundary condition to folder 'boundaryData'.")
        args = [time_step, boundary_inlet_path, alpha_name, n_faces, U_inlet, VOF_inlet, time_inlet]
        partial_write = partial(self.__write_time_step, args=args)
        n_workers = max(8, os.cpu_count())
        self.logger.info(f"Using {n_workers} threads for writing routine.")
        time_steps = [i for i in range(len(time_inlet))]
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            executor.map(partial_write, time_steps)
        self.logger.info("Boundary condition was successfully saved in 'boundaryData'.")


    # === Private functions ===
    def __write_header(self, file_loc: str, class_name: str, object_name: str):
        with open(file_loc, 'w') as f:
            f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
            f.write(r'| =========                 |                                                 |'+"\n")
            f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
            f.write(r'|  \\    /   O peration     | Version:  x                                   |'+"\n")
            f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
            f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
            f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
            f.write(r'FoamFile'+"\n")
            f.write(r'{'+"\n")
            # f.write('\t version \t\t x;'+"\n")
            f.write('\t format \t\t ascii;'+"\n")
            f.write('\t class \t\t ' + class_name + ';'+"\n")
            f.write('\t object \t\t ' + object_name + ';'+"\n");
            f.write('}'+"\n")
            f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
            f.write("\n")

    def __write_footer(self, file_loc: str):
        with open(file_loc, 'a+') as f:
            f.write("\n")
            f.write(r'// ************************************************************************* //' + "\n")

    def __write_time_step(self, i, args):
        # arguments
        time_step = args[0]
        boundary_inlet_path = args[1]
        alpha_name = args[2]
        n_faces = args[3]
        U_inlet = args[4]
        VOF_inlet = args[5]
        time_inlet = args[6]


        time_i = time_inlet[i]
        trunc_time_i = truncate(time_i, time_step)

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

        # Write OpenFOAM-header
        self.__write_header(U_out_path, "vectorAverageField", "values")
        self.__write_header(VOFw_out_path, "scalarAverageField", "values")

        # write content
        with open(U_out_path, 'a+') as f:
            f.write('//Data points'+'\n')
            f.write(str(n_faces)+'\n'+'('+'\n')
            for j in range(n_faces):
                f.write('(' + str(U_inlet[j, i, 0]) + ' ' + str(U_inlet[j, i, 1]) + ' ' + str(U_inlet[j, i, 2]) + ') \n')
            f.write(')'+'\n')

        with open(VOFw_out_path, 'a+') as f:
            f.write('//Data points'+'\n')
            f.write(str(n_faces)+'\n'+'('+'\n')
            for j in range(n_faces):
                f.write(str(VOF_inlet[j,i,0])+'\n')
            f.write(')'+'\n')

        # Write OpenFOAM-footer
        self.__write_footer(U_out_path)
        self.__write_footer(VOFw_out_path)