from utils import os, yaml

# TODO additional functions for openfoam foundation 

def create_write_cell_centres_file(file_path):
    f = open(file_path, 'w')

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
    
    f.close()

def create_write_cell_areas_file(file_path):
    f = open(file_path, 'w')

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

    # code execution function
    f.write("\tcodeExecute\n")
    f.write("\t#{\n")
    # code content
    f.write("\t\tInfo << \"Execute writeFaceAreas\" << endl;\n")
    f.write("\t\tlabel patchId = mesh().boundaryMesh().findPatchID(\"inlet\");\n")
    f.write("\t\tsurfaceScalarField faceAreas\n")
    f.write("\t\t(\n")

    f.write("\t\t\tIOobject\n")
    f.write("\t\t\t(\n")

    f.write("\t\t\t\t\"area\",\n")
    f.write("\t\t\t\tmesh().time().timeName(),\n")
    f.write("\t\t\t\tmesh(),\n")
    f.write("\t\t\t\tIOobject::NO_READ,\n")
    f.write("\t\t\t\tIOobject::NO_WRITE,\n")
    f.write("\t\t\t\tIOobject::NO_REGISTER\n")

    f.write("\t\t\t),\n")
    f.write("\t\t\tmesh().magSf()\n")

    f.write("\t\t);\n")
    f.write("\t\tfaceAreas.write();\n")

    f.write("\t#};\n")
    f.write("}\n")
     
    f.close()

def modify_control_dict(file_path):
    functions_line_exist = [-1, False]
    with open(file_path, "r") as f:
        buf = f.readlines()

        # check if functions exist
        for line in buf:
            functions_line_exist[0] += 1
            if line=="functions\n":
                functions_line_exist[1] = True
                break
    
    # write include FOs
    if not functions_line_exist[1]:

        f = open(file_path, "r").readlines()
        f[-1] = "functions\n"
        f.append("{\n")
        f.append("\t#include \"FOs/FOwriteCellCentres\"\n")
        f.append("\t#include \"FOs/FOwriteCellAreas\"\n")
        f.append("}\n")
        f.append("// ************************************************************************* //")

        open(file_path, "w").writelines(f)
    else:
        f = open(file_path, "r").readlines()

        # find functions line
        functions_idx = f.index("functions\n")
        functions_end_idx = f[functions_idx:].index("}\n")
        
        # check if includes already added
        fo_writecellcentres = "\t#include \"FOs/FOwriteCellCentres\"\n"
        fo_writecellareas = "\t#include \"FOs/FOwriteCellAreas\"\n"
        fo_writecellcentres_exist = fo_writecellcentres in f[functions_idx:]
        fo_writecellareas_exist = fo_writecellareas in f[functions_idx:]
            

        temp = f[:functions_idx+2] 
        if not fo_writecellcentres_exist:
            temp.append(fo_writecellcentres)
        if not fo_writecellareas_exist:
            temp.append(fo_writecellareas)
        temp = temp + f[functions_idx+2:]

        open(file_path, "w").writelines(temp)
                

def create_postprocess_functions(config):
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    case_path = config["case_path"]

    # create FO directory and files if not exist
    fo_dir_name = "FOs"
    fo_dir_path = os.path.join(case_path, "system", fo_dir_name)
    if not os.path.exists(fo_dir_path):
        os.mkdir(fo_dir_path)

    # write cell centres
    writecellcentres_path = os.path.join(fo_dir_path, "FOwriteCellCentres")
    if not os.path.exists(writecellcentres_path):
        create_write_cell_centres_file(writecellcentres_path)

    # write cell areas
    writecellareas_path = os.path.join(fo_dir_path, "FOwriteCellAreas")
    if not os.path.exists(writecellareas_path):
        create_write_cell_areas_file(writecellareas_path)
    
    # include FOs to controlDict
    controldict_path = os.path.join(case_path, "system", "controlDict")
    if os.path.exists(controldict_path):
        modify_control_dict(controldict_path)
    else:
        raise RuntimeError("controlDict does not exist.")
