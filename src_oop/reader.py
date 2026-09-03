from utils import np, os, subprocess, linecache

class Reader():
    TOLERANCE = 0.01

    def __init__(self, config, data):
        cfd_prog = str(config["packages"]["cfd_program"])
        cfd_version = str(config["packages"]["cfd_version"])
        self.cfd_module = f"{cfd_prog}/{cfd_version}"
        self.dimension = int(config["cfd"]["dimension"])
        self.inlet_name = str(config["cfd"]["inlet_name"])
        self.time_start = str(config["model"]["time"]["start"])
        self.data = data

        self.output_path = str(config.get("_output_path"))
        self.openfoam_type = str(config.get("_openfoam_type"))

        # protected variables
        self.inlet_cell_data_names = None
        if self.openfoam_type=="org":
            self.inlet_cell_data_names = ["Ccx", "Ccy", "Ccz", "area"]
        else:
            self.inlet_cell_data_names = ["Cx", "Cy", "Cz", "area"]

    def initialize(self):
        '''initialize data storage for inlet geometry'''
        self.inlet_faces = None
        self.inlet_normal = None

    def run(self):
        self._run_openfoam_postprocess()
        self._read_cell_variables()
        self._check_face_list()

        if self.dimension==3:
            self._calculate_inlet_normal_3d()
        elif self.dimension==2:
            self._calculate_inlet_normal_2d()
        else:
            raise RuntimeError(f"Number of dimensions {self.dimension} is invalid. It must either be 2 or 3.")

        self._pass_data()
        self._save_inlet_faces()

    def get_data(self):
        return self.inlet_faces, self.inlet_normal

    # ===== Protected functions =====
    def _run_openfoam_postprocess(self):
        openfoam_type = self.openfoam_type
        cfd_module = self.cfd_module
        time_start = self.time_start

        if openfoam_type=="com":
            command = "ml " + cfd_module + "; source $FOAM_BASH; postProcess -time " + time_start + ";"
        if openfoam_type=="org":
            command = "ml " + cfd_module + "; source $FOAM_BASH; foamPostProcess; foamPostProcess -func writeCellCentres;"

        subprocess.run(command, cwd=os.getcwd(), shell=True)

    def _read_cell_variables(self):
        # attributes
        cwd = os.getcwd()
        inlet_name = self.inlet_name
        time_start = self.time_start

        # Loading inlet cell coordinates and face areas into Python
        for i, file_i in enumerate(self.inlet_cell_data_names):
            source_file = f"{cwd}/{time_start}/{file_i}"
            if not os.path.exists(source_file):
                raise RuntimeError(f"{source_file} not found")

            try:
                command = "grep -nr " + inlet_name + " " + source_file + " | cut -d : -f 1 > lineNr"
                subprocess.run(command, cwd=cwd, shell=True)
                lineNameNr = int(open(cwd + "/lineNr", 'r').readline())
                lineStartNr = lineNameNr + 6  # In case of non-uniform list, this is where the list of values in the source_file starts
                rowsNrIndex = lineNameNr + 4  # On this line, the number of cell centers on the inlet is stated

                command = "awk NR==" + str(rowsNrIndex) + " " + source_file + " > rowsNr"
                subprocess.run(command, cwd=cwd, shell=True)
                rowsNr = int(open(cwd + "/rowsNr", 'r').readline())
                os.remove(os.path.join(cwd, "lineNr"))
                os.remove(os.path.join(cwd, "rowsNr"))
                tempCoordFile = np.ones([rowsNr, 1]) * float("inf")
                for j in np.arange(rowsNr):
                    tempCoordFile[j, 0] = float(linecache.getline(source_file, lineStartNr + j))
            # If the source-file has uniform coordinate
            except ValueError:
                # If 'rowsNr' not local variable
                if 'rowsNr' not in locals():
                    center_y_file = cwd + "/" + time_start + self.inlet_cell_data_names[1]

                    command = "grep -nr " + inlet_name + " " + center_y_file + " | cut -d : -f 1 > lineNr"
                    subprocess.run(command, cwd=cwd, shell=True)
                    lineNameNr_CF = int(open(cwd+"/lineNr",'r').readline())
                    rowsNrIndex_CF = lineNameNr_CF+4 # On this line, the number of cell centers on the inlet is stated

                    command = "awk NR==" + str(rowsNrIndex_CF) + " " + center_y_file + " > rowsNr"
                    subprocess.run(command, cwd=cwd, shell=True)
                    rowsNr = int(open(cwd+"/rowsNr", 'r').readline())
                    os.remove(os.path.join(cwd, "lineNr"))
                    os.remove(os.path.join(cwd, "rowsNr"))

                indexUV = lineNameNr + 3
                command = "awk NR==" + str(indexUV) + " " + source_file + " > unifValue"
                subprocess.run(command, cwd=cwd, shell=True)
                unifValue = float(open(cwd+"/unifValue", 'r').readline().split()[-1][0:-1])
                os.remove(os.path.join(cwd, "unifValue"))
                tempCoordFile = np.ones([rowsNr, 1])*float("inf")
                for j in np.arange(rowsNr):
                    tempCoordFile[j, 0] = unifValue

            if i == 0:
                self.inlet_faces = np.ones([rowsNr, 5]) * float("inf")  # ID - X - Y - Z - Area
                self.inlet_faces[:, 0] = np.arange(rowsNr)

            self.inlet_faces[:, (i + 1)] = tempCoordFile[:, 0]

    def _check_face_list(self):
        # Check that all values are inserted correctly
        is_facelist_too_large = self.inlet_faces > 1e16
        is_facelist_wrong = is_facelist_too_large.nonzero()[0]
        if is_facelist_wrong:
            raise RuntimeError("Not all values are correctly read.")

    def _calculate_inlet_normal_3d(self):
        cwd = os.getcwd()

        point1 = self.inlet_faces[0, 1:4]
        point2 = self.inlet_faces[1, 1:4]
        i = 2
        point3 = self.inlet_faces[i, 1:4]
        enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
        denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))
        while (enum / denom) < self.TOLERANCE:
            i = i + 1
            point3 = self.inlet_faces[i, 1:4]
            enum = np.linalg.norm(np.cross(point2 - point1, point3 - point1))
            denom = (np.linalg.norm(point2 - point1) * np.linalg.norm(point3 - point1))

        normal_vec = np.cross(point2 - point1, point3 - point1)
        normal_inlet = normal_vec / np.linalg.norm(normal_vec)

        # Need one more point from the domain to determine the correct orientation of the inlet normal
        command = "grep -nr '(' constant/polyMesh/points | head -n 1 | cut -d : -f 1 > lineNr"
        subprocess.run(command, cwd=cwd, shell=True)
        lineNameNr = int(open(cwd+"/lineNr", 'r').readline())
        lineNr = lineNameNr+1  # First point that is defined
        os.remove(os.path.join(cwd, "lineNr"))
        i = 0
        points_path = os.path.join(cwd, "constant", "polyMesh", "points")
        with open(points_path, 'r') as f:
            while i < lineNr:
                i = i+1
                f.readline()
            point_inside_domain = np.double(np.array(f.readline()[1:-2].split(" ")))

            #find point_inside_domain which is not in the inlet plane
            while abs(np.dot(point_inside_domain-point1, normal_inlet))/(np.linalg.norm(point_inside_domain-point1)) < self.TOLERANCE:
                point_inside_domain = np.double(np.array(f.readline()[1:-2].split(" ")))

            # switch sign if normal_inlet pointing outside
            if np.dot(point_inside_domain-point1, normal_inlet) < 0:
                normal_inlet = (-1.0)*normal_inlet

        self.inlet_normal = normal_inlet

    def _calculate_inlet_normal_2d(self):
        normal_inlet = np.zeros([3])
        print("The normal to the inlet cannot be calculated directly. Please input the x-, y- and z-coordinates of the normal vector.")
        axis = ["x", "y", "z"]
        for i, ax_i in enumerate(axis):
            is_normal_correct = False
            while not is_normal_correct:
                try:
                    print(f"Please provide the {ax_i}-component of the normal vector: ")
                    temp_normal = float(input())
                    is_normal_correct = True
                except ValueError:
                    print("Input values is not expected type. Please provide a float.")
            normal_inlet[i] = temp_normal
        normal_inlet = (1/np.linalg.norm(normal_inlet))*normal_inlet

        self.inlet_normal = normal_inlet

    def _pass_data(self):
        self.data["inlet_faces"] = self.inlet_faces
        self.data["inlet_normal"] = self.inlet_normal

    def _save_inlet_faces(self):
        # Save inlet and normal in Python Numpy-array format
        np.save(os.path.join(self.output_path, "inlet_faces.npy"), self.inlet_faces)
        np.save(os.path.join(self.output_path, "inlet_normal.npy"), self.inlet_normal)