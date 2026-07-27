from utils import os, np, plt, re, pd
from openfoam_io import get_int, get_vector_array

# regex patterns
cofr_pattern = r"# CofR.*"
float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
int_pattern = r'[+-]?\d+'
delimiter = r'[\s\n]+'

## USER INPUT
variable_map = {
    "force": "Force [N]",
    "moment": "Moment [Nm]",
    "U": "Velocity [m/s]"
    }

def read_cofr_location(probes_path):
    files = os.listdir(probes_path)

    cofr_loc = None
    # read first probe data file to get probe location
    for file_i in files[0]:
        fi_path = os.path.join(probes_path, file_i)
        for line in open(fi_path, "r"):
            if re.search(cofr_pattern, line):
                cofr_loc = get_vector_array(line)[0]
            else:
                break

    return cofr_loc

def read_data(probes_path) -> dict:
    files = os.listdir(probes_path)
    colnames = ["time", "total_x", "total_y", "total_z",
                "pressure_x", "pressure_y", "pressure_z",
                "viscous_x", "viscous_y", "viscous_z",
                "NaN" # required, dropped later
                ]

    force_data = {}
    for file_i in files:
        fi_path = os.path.join(probes_path, file_i)

        delim = r"\s+"
        data = pd.read_csv(fi_path,
                           header=3, # number of header text in both forces.dat and moment.dat
                           delimiter=delim,
                           names=colnames,
                           )

        # drop column which contain only NaNs
        data = data.drop("NaN", axis=1)
        force_data[file_i[:-4]] = data

    return force_data

def plot_data(data: dict, var_name: str, key_type=None, plot_range=None, marker=None):
    '''
    Plot data[var_name].

    :param data:        dictionary containing data.
    :param var_name:    variable name. must be in data.keys()
    :param key_type:    list of keyword to filter columns. Use None to plot all columns.
    '''
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, layout='constrained')

    # time variables
    time = data[var_name]["time"]
    dt = time[1] - time[0] # assume const timestep

    force_keys = data[var_name].keys()[1:]
    target_keys = force_keys
    if key_type:
        target_keys = []
        for key in force_keys:
            for match in key_type:
                if match in key:
                    target_keys.append(key)

    target_keys = list(set(target_keys))
    target_keys.sort()

    for key in target_keys:
        var = data[var_name][key]

        if plot_range:
            x = time[plot_range]
            y = var[plot_range]
        else:
            x = time
            y = var

        # time space plot
        ax0.plot(x, y, label=key, marker=marker)
        ax0.set_ylabel(variable_map[var_name])
        ax0.set_xlabel("Time [s]")
        ax0.legend()
        ax0.grid(True)

        # PSD
        ax1.psd(y-np.mean(y), NFFT=2**10, Fs=1/dt, label=key, scale_by_freq=True)
        ax1.set_ylabel("Power Spectral Density [dB/Hz]")
        ax1.set_xlabel("Frequency [Hz]")
        # ax1.legend()
        ax1.grid(True)

        # Spcetral analysis
        ax2.magnitude_spectrum(y-np.mean(y), Fs=1/dt)
        ax2.set_ylabel("Magnitude spectrum")
        ax2.set_xlabel("Frequency [Hz]")
        # ax2.legend()
        ax2.grid(True)

if __name__=="__main__":
    # get paths
    case_path = os.getcwd()
    probes_path = os.path.join(case_path, "postProcessing", "forces", "0")

    force_data = read_data(probes_path)
    time = force_data["force"]["time"]
    dt = time[1] - time[0] # assume const timestep
    plot_range = None #slice(int(0.2/dt), int(0.3/dt))

    plot_data(force_data, var_name="force", key_type=["total"],
              plot_range=plot_range, marker=None)
    plot_data(force_data, var_name="moment", key_type=["total"],
              plot_range=plot_range, marker=None)
    plt.show()

    print("plot_forces script ended successfully.")
