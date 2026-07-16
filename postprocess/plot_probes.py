from utils import os, plt, re, pd
from openfoam_io import get_int, get_vector_array

# regex patterns
probe_loc_pattern = r"# Probe.*"
float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
int_pattern = r'[+-]?\d+'
delimiter = r'[\s\n]+'

## USER INPUT
variable_map = {
    "p": "Pressure [Pa]",
    "alpha.Water": "Alpha water [-]",
    "U": "Velocity [m/s]"
    }

def read_probe_location(probes_path):
    files = os.listdir(probes_path)

    probe_loc = {}
    # read first probe data file to get probe location
    for file_i in files[0]:
        fi_path = os.path.join(probes_path, file_i)
        for line in open(fi_path, "r"):
            if re.search(probe_loc_pattern, line):
                probe_i = get_int(line, "# Probe")
                probe_loc[f"P{probe_i}"] = get_vector_array(line)[0]
            else:
                break

    return probe_loc

def read_probe_data(probe_loc, probes_path):
    files = os.listdir(probes_path)
    n_probes = len(probe_loc.keys())
    colnames = ["time"] + [f"P{i}" for i in range(n_probes)]

    probe_data = {}
    for file_i in files:
        fi_path = os.path.join(probes_path, file_i)

        delim = r"\s+"
        if file_i=="U":
            delim = r"\s+(?=\()"

        data = pd.read_csv(fi_path,
                           header=n_probes+1,
                           delimiter=delim,
                           names=colnames,
                           )
        probe_data[file_i] = data

    return probe_data

def plot_probe_data(probe_data, var_name, probe_idx=None):
    fig, ax = plt.subplots()
    time = probe_data[var_name]["time"]
    if not probe_idx:
        n_probes = range(len(probe_data[var_name].keys())-1)
    else:
        n_probes = probe_idx
    for i in n_probes:
        probe_i = f"P{i}"
        var = probe_data[var_name][probe_i]
        ax.plot(time, var, label=probe_i, marker="o")
        ax.set_ylabel(variable_map[var_name])
        ax.set_xlabel("Time [s]")
        ax.legend()

if __name__=="__main__":
    # get paths
    case_path = os.getcwd()
    probes_path = os.path.join(case_path, "postProcessing", "probes", "0")

    probe_loc = read_probe_location(probes_path)
    probe_data = read_probe_data(probe_loc, probes_path)

    # test print
    # temp = probe_data["p"]["time"].to_numpy()
    # temp2 = probe_data["U"]["P1"].to_numpy()

    plot_probe_data(probe_data, var_name="p", probe_idx = None)
    plot_probe_data(probe_data, var_name="alpha.Water", probe_idx = [0,2,4])
    plt.show()

    print("0")
