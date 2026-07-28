from utils import os, plt, re, pd, np
from openfoam_io import get_int, get_vector_array

# regex patterns
probe_loc_pattern = r"# Probe.*"
delimiter = r'[\s\n]+'

## USER INPUT
PSD_XRANGE = [0, 150] # in Hz
PSD_WINDOW = [0.05, 0.1] # in s
VARIABLE_MAP = {
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

def plot_probe_data(probe_loc, probe_data, var_name,
                    probe_idx=None, scale=1000, psd_xrange=[None, None], psd_window=None,
                    ):
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, layout='constrained')

    # time variables
    time = probe_data[var_name]["time"]
    dt = time[1] - time[0] # assume const timestep

    if not probe_idx:
        n_probes = range(len(probe_data[var_name].keys())-1)
    else:
        n_probes = probe_idx

    for i in n_probes:
        # define plot labels
        probe_i = f"P{i}"
        probe_loc_i = probe_loc[probe_i] * scale
        temp = np.array2string(probe_loc_i, precision=2, separator=",", suppress_small=True)
        label = probe_i + " " + temp

        # plot variables
        var = probe_data[var_name][probe_i]
        if psd_window:
            x = time[psd_window]
            y = var[psd_window]
        else:
            x = time
            y = var

        marker = "o"
        if len(var)>100:
            marker = None
        ax0.plot(x, y, label=label, marker=marker)
        ax0.set_ylabel(VARIABLE_MAP[var_name])
        ax0.set_xlabel("Time [s]")
        ax0.grid(True)
        # ax0.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        # ax0.legend(ncol=4, bbox_to_anchor=(0.5,-0.5), loc='lower center', edgecolor='w')

        # PSD
        interp_points = 1000 #len(y)
        print(f"PSD using {interp_points} interpolation points")
        ax1.psd(y-np.mean(y), NFFT=interp_points, Fs=1/dt, scale_by_freq=True)
        ax1.set_ylabel("Power Spectral Density [dB/Hz]")
        ax1.set_xlabel("Frequency [Hz]")
        ax1.set_xlim(psd_xrange[0], psd_xrange[1])
        # ax1.legend()

        # Spectral analysis
        ax2.magnitude_spectrum(y-np.mean(y), Fs=1/dt)
        ax2.set_ylabel("Magnitude spectrum")
        ax2.set_xlabel("Frequency [Hz]")
        ax2.set_xlim(psd_xrange[0], psd_xrange[1])
        # ax2.legend()

    fig.legend(ncol=4, loc='outside upper right')

if __name__=="__main__":
    # get paths
    case_path = os.getcwd()
    probes_path = os.path.join(case_path, "postProcessing", "probes", "0")

    probe_loc = read_probe_location(probes_path)
    probe_data = read_probe_data(probe_loc, probes_path)

    time = probe_data["p"]["time"]
    dt = time[1] - time[0] # assume const timestep
    psd_window = slice(int(PSD_WINDOW[0]/dt), int(PSD_WINDOW[1]/dt))

    plot_probe_data(probe_loc, probe_data,
                    var_name="p", probe_idx = [0,1,2],
                    psd_xrange=PSD_XRANGE, psd_window=psd_window)
    plot_probe_data(probe_loc, probe_data,
                    var_name="alpha.Water", probe_idx = [0,2,4],
                    psd_xrange=PSD_XRANGE, psd_window=psd_window)
    plt.show()
    print("0")
