from utils import os, re

def create_write_cell_centres(file_path):
    with open(file_path, "w") as f:
        f.write("writeCellCentres1\n")
        f.write("{\n")

        # mandatory entries
        f.write("\ttype writeCellCentres;\n")
        f.write("\tlibs (fieldFunctionObjects);\n")

        f.write("\n")

        # optional entries
        f.write("\tenabled true;\n")
        f.write("\tlog true;\n")
        f.write("\ttimeStart 0;\n")
        f.write("\ttimeEnd 0;\n")
        f.write("\texecuteControl writeTime;\n")
        f.write("\texecuteInterval -1;\n")
        f.write("\twriteControl writeTime;\n")
        f.write("\twriteInterval -1;\n")

        f.write("}\n")

def create_write_cell_areas(file_path, openfoam_type):
    with open(file_path, 'w') as f:
        f.write("writeCellAreas1\n")
        f.write("{\n")

        f.write("\ttype coded;\n")
        f.write("\tlibs (\"libutilityFunctionObjects.so\");\n")
        f.write("\tname writeCellAreas;\n")

        f.write("\n")

        # optional entries
        f.write("\tenabled true;\n")
        f.write("\tlog true;\n")
        f.write("\ttimeStart 0;\n")
        f.write("\ttimeEnd 0;\n")
        f.write("\texecuteControl writeTime;\n")
        f.write("\texecuteInterval -1;\n")
        f.write("\twriteControl writeTime;\n")
        f.write("\twriteInterval -1;\n")

        f.write("\n")

        # write code options necessary for org openfoam
        if openfoam_type == "org":
            f.write("\tcodeInclude\n")
            f.write("\t#{\n")
            f.write("\t\t#include \"volFields.H\"\n")
            f.write("\t\t#include \"surfaceFields.H\"\n")
            f.write("\t#};\n")
            f.write("\n")

        # code execution function
        f.write("\tcodeExecute\n")
        f.write("\t#{\n")
        # code content
        f.write("\t\tInfo << \"Execute writeFaceAreas\" << endl;\n")
        # f.write("\t\tlabel patchId = mesh().boundaryMesh().findPatchID(\"inlet\");\n")
        f.write("\t\tsurfaceScalarField faceAreas\n")
        f.write("\t\t(\n")

        f.write("\t\t\tIOobject\n")
        f.write("\t\t\t(\n")

        f.write("\t\t\t\t\"area\",\n")
        if openfoam_type == "com":
            f.write("\t\t\t\tmesh().time().timeName(),\n")
        if openfoam_type == "org":
            f.write("\t\t\t\tmesh().time().name(),\n")
        f.write("\t\t\t\tmesh(),\n")
        f.write("\t\t\t\tIOobject::NO_READ,\n")

        if openfoam_type == "com":
            f.write("\t\t\t\tIOobject::NO_WRITE,\n")
            f.write("\t\t\t\tIOobject::NO_REGISTER\n")
        if openfoam_type == "org":
            f.write("\t\t\t\tIOobject::AUTO_WRITE\n")

        f.write("\t\t\t),\n")
        f.write("\t\t\tmesh().magSf()\n")

        f.write("\t\t);\n")
        f.write("\n")

        f.write("\t\tInfo << \"\tWriting cell face areas [area] to 0\" << endl;\n")
        f.write("\t\tfaceAreas.write();\n")

        f.write("\t#};\n")
        f.write("}\n")

class Preprocessor():
    def __init__(self, config):
        self.cfd_program = str(config["packages"]["cfd_program"])
        self.cfd_version = str(config["packages"]["cfd_version"])

        self.openfoam_type = str(config.get("_openfoam_type"))

        # paths
        self.cwd = os.getcwd()
        self.fo_dir_path = None

    def run(self):
        # run routine
        self._create_FO_directory()
        self._write_cell_variables()
        self._modify_controlDict()

    # ===== Protected functions =====
    def _create_FO_directory(self):
        fo_dir_name = "FOs"
        self.fo_dir_path = os.path.join(self.cwd, "system", fo_dir_name)
        if not os.path.exists(self.fo_dir_path):
            os.mkdir(self.fo_dir_path)

    def _write_cell_variables(self):
        writecellcentres_path = os.path.join(self.fo_dir_path, "FOwriteCellCentres")
        if not os.path.exists(writecellcentres_path):
            create_write_cell_centres(writecellcentres_path)

        writecellareas_path = os.path.join(self.fo_dir_path, "FOwriteCellAreas")
        if not os.path.exists(writecellareas_path):
            create_write_cell_areas(writecellareas_path, self.openfoam_type)

    def _modify_controlDict(self):
        '''Include FOs to controlDict.'''
        openfoam_type = self.openfoam_type
        controldict_path = os.path.join(self.cwd, "system", "controlDict")
        if not os.path.exists(controldict_path):
            raise RuntimeError(f"controlDict does not exist in {controldict_path}.")

        functions_line_exist = [-1, False]
        with open(controldict_path, "r") as f:
            buf = f.readlines()

            # check if functions exist
            for line in buf:
                functions_line_exist[0] += 1
                if line=="functions\n":
                    functions_line_exist[1] = True
                    break

        # write include FOs
        if not functions_line_exist[1]:

            f = open(controldict_path, "r").readlines()
            f[-1] = "functions\n"
            f.append("{\n")
            if openfoam_type=="com":
                f.append("\t#include \"FOs/FOwriteCellCentres\"\n")
            f.append("\t#include \"FOs/FOwriteCellAreas\"\n")
            f.append("}\n")
            f.append("// ************************************************************************* //")

            open(controldict_path, "w").writelines(f)
        else:
            f = open(controldict_path, "r").readlines()

            # find functions line
            functions_idx = f.index("functions\n")
            # functions_end_idx = f[functions_idx:].index("}\n")

            # check if includes already added
            fo_writecellcentres = "\t#include \"FOs/FOwriteCellCentres\"\n"
            fo_writecellareas = "\t#include \"FOs/FOwriteCellAreas\"\n"
            fo_writecellcentres_exist = fo_writecellcentres in f[functions_idx:]
            fo_writecellareas_exist = fo_writecellareas in f[functions_idx:]


            temp = f[:functions_idx+2]
            if (not fo_writecellcentres_exist) & (openfoam_type=="com"):
                temp.append(fo_writecellcentres)
            if not fo_writecellareas_exist:
                temp.append(fo_writecellareas)
            temp = temp + f[functions_idx+2:]

            open(controldict_path, "w").writelines(temp)