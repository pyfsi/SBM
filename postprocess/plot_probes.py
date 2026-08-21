from utils import os, plt, re, pd, np, sps
from openfoam_io import get_int, get_vector_array

# regex patterns
probe_loc_pattern = r"# Probe.*"
delimiter = r'[\s\n]+'

## USER INPUT
NFFT = 1 << 11          # points for Welch algorithm
F_SAMPLING = None        # sampling frequency for Welch algorithm
PSD_XRANGE = [0, 500]   # PSD plot x-axis range in Hz
PSD_WINDOW = [0.0, 0.05] # window for PSD function. Ignore all data points outside this range.
PROBE_IDX = None  # probe indices

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
    fig, (ax0, ax1) = plt.subplots(2, 1, layout="constrained")
    fig.tight_layout(h_pad=4, rect=[0.02, 0.05, 0.83, 0.98])

    # time variables
    time = probe_data[var_name]["time"]
    dt = time[1] - time[0] # assume const timestep
    fs = 1/dt
    if F_SAMPLING:
        fs = F_SAMPLING

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

        # convert to numpy
        x = x.to_numpy()
        y = y.to_numpy()

        marker = "o"
        ax0.plot(x, y, label=label, marker=marker)
        ax0.set_ylabel(VARIABLE_MAP[var_name])
        ax0.set_xlabel("Time [s]")
        ax0.grid(True)
        ax0.ticklabel_format(axis='y', style='sci', scilimits=(0,5))

        # PSD
        freq, pxx = sps.welch(y, fs=fs, nfft=NFFT,)
        ax1.plot(freq, pxx, marker=marker)
        ax1.set_ylabel("PSD [m*m/Hz]")
        ax1.set_xlabel("Frequency [Hz]")
        ax1.set_xlim(psd_xrange[0], psd_xrange[1])
        ax1.grid(True)
        ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,5))

        # legend
        fig.legend(ncol=1,
                   bbox_to_anchor=(1.0, 1.0),
                   loc='upper right',
                   edgecolor='k',
                   title="Probe locations in mm",
                   fontsize='small',
                   fancybox=True)

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
                    var_name="p", probe_idx = PROBE_IDX,
                    psd_xrange=PSD_XRANGE, psd_window=psd_window)
    plot_probe_data(probe_loc, probe_data,
                    var_name="alpha.Water", probe_idx = PROBE_IDX,
                    psd_xrange=PSD_XRANGE, psd_window=psd_window)
    plt.show()
