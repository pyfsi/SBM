from utils import np, sys, os, random, logging, plt
from utils import PI, SBM_OUTPUT


class Plotter():
    def __init__(self, config):
        # store internal variables
        self.mg_min = float(config["sbm"]["mass_bubble_min"])
        self.mg_max = float(config["sbm"]["mass_bubble_max"])

        # storage for bubble parameter from sampling alogrithm
        self.sample_data_idx = 0
        self.sample_data = {}
        self.sample_data["face_idx"] = [None] * 10000
        self.sample_data["time_idx"] = [None] * 10000

        # distribution parametrs (bins per bubble variable, counts)
        n_bins = 100
        mass_bin_width = (self.mg_max - self.mg_min) / n_bins
        self.distribution = {}
        self.distribution["mass"] = (
            np.arange(self.mg_min/mass_bin_width, self.mg_max/mass_bin_width+1, 1) * mass_bin_width,
            np.zeros(n_bins)
        )
        self.distribution["mass_sample"] = (
            np.arange(self.mg_min/mass_bin_width, self.mg_max/mass_bin_width+1, 1) * mass_bin_width,
            np.zeros(n_bins)
        )

        # relative absolute difference
        self.residual_idx = 0
        self.residual = [None] * 10000

        # fig and axes
        self.fig, self.axs = plt.subplots(2, 2, layout='constrained')
        # ax 0,0 settings
        self.axs[0,0].set_ylabel('Residual [-]')
        self.axs[0,0].set_xlabel('Number of iterations [-]')
        self.axs[0,0].set_ylim(0.0, 1.0)
        # ax 1,0 settings
        self.axs[1,0].set_ylabel('Count [-]')
        self.axs[1,0].set_xlabel('Bubble mass [kg]')
        # ax 0,1 settings
        self.axs[0,1].set_ylabel('Face index [-]')
        self.axs[0,1].set_xlabel('Time index [-]')
        # ax 1,1 settings
        self.axs[1,1].set_ylabel('Count [-]')
        self.axs[1,1].set_xlabel('Sample mass [kg]')

    def _clear_figures(self):
        # clear lines and patches
        for ax_i in self.axs[:,:]:
            for ax_ij in ax_i[:]:
                for line in ax_ij.lines:
                    line.remove()
                for patch in ax_ij.patches:
                    patch.remove()
                for text in ax_ij.texts:
                    text.remove()

    def update_data(self, face_idx, time_idx):
        # store sample data
        data_idx = self.sample_data_idx
        self.sample_data["face_idx"][data_idx] = face_idx
        self.sample_data["time_idx"][data_idx] = time_idx

        # increment index
        self.sample_data_idx += 1

        # resize sample data array if it is half full
        n_curr_sample_data = len(self.sample_data["face_idx"])
        if self.sample_data["face_idx"].index(None) > (n_curr_sample_data//2):
            self.sample_data["face_idx"] += [None] * n_curr_sample_data
            self.sample_data["time_idx"] += [None] * n_curr_sample_data

    def update_sample_distribution(self, mass_sample):
        # get bin index for distribution plot
        mass_sample_idx = np.digitize([mass_sample], self.distribution["mass_sample"][0]) - 1

        # increment count at bin
        self.distribution["mass_sample"][1][mass_sample_idx] += 1

    def update_distribution(self, mass_bubble_cells):
        # get bin index for distribution plot
        mass_idx = np.digitize([mass_bubble_cells], self.distribution["mass"][0]) - 1

        # increment count at bin
        self.distribution["mass"][1][mass_idx] += 1

    def update_residual(self, mass_per_block, mass_inserted):
        # define residual as relative absolute difference
        residual = abs(mass_per_block - mass_inserted) / mass_per_block
        self.residual[self.residual_idx] = residual

        # move array pointer idx
        self.residual_idx += 1

        # resize residual array if it is half full
        n_curr_residual = len(self.residual)
        if self.residual.index(None) > (n_curr_residual//2):
            self.residual += [None] * n_curr_residual

    def _plot_hist(self,):
        # plot settings template
        box_props = dict(boxstyle='round', edgecolor='black', facecolor="white", alpha=0.9)
        bar_kwargs = {"x":None, "height":None,
                      "edgecolor":"black", "facecolor":None, "width":None,}
        txt_kwargs = {"x":0.02, "y":0.93, "s":"",
                      "transform":None, "fontsize":"15","bbox":box_props}

        # === plot defined bubble mass distribution ===
        # calculate distribution variables
        n_bubbles = np.sum(self.distribution["mass"][1])
        mass_bar_x = self.distribution["mass"][0][:-1] + 0.5 * np.diff(self.distribution["mass"][0][:])

        # assign to kwargs
        bar_kwargs["x"] = mass_bar_x
        bar_kwargs["height"] = self.distribution["mass"][1]
        bar_kwargs["facecolor"] = "green"
        bar_kwargs["width"] = self.distribution["mass"][0][1]-self.distribution["mass"][0][0]
        txt_kwargs["transform"] = self.axs[1,0].transAxes
        txt_kwargs["s"] = f"n={int(n_bubbles)}"

        # plot
        self.axs[1,0].bar(**bar_kwargs)
        self.axs[1,0].text(**txt_kwargs)

        # === plot sample mass distribution ===
        # calculate distribution variables
        n_samples = np.sum(self.distribution["mass_sample"][1])
        sample_mass_bar_x = self.distribution["mass_sample"][0][:-1] + 0.5 * np.diff(self.distribution["mass_sample"][0][:])

        # assign to kwargs
        bar_kwargs["x"] = sample_mass_bar_x
        bar_kwargs["height"] = self.distribution["mass_sample"][1]
        bar_kwargs["facecolor"] = "red"
        bar_kwargs["width"] = self.distribution["mass_sample"][0][1]-self.distribution["mass_sample"][0][0]
        txt_kwargs["transform"] = self.axs[1,1].transAxes
        txt_kwargs["s"] = f"n={int(n_samples)}"

        # plot
        self.axs[1,1].bar(**bar_kwargs)
        self.axs[1,1].text(**txt_kwargs)

    def plot_bubble_insertion(self):
        # remove previous figures
        self._clear_figures()

        # plot residual
        x_ticks = np.arange(1, len(self.residual)+1, 1)
        self.axs[0,0].plot(
            x_ticks,
            np.array(self.residual),
            color="black",
            marker="o",
        )

        # plot face and time index
        self.axs[0,1].scatter(
            self.sample_data["time_idx"],
            self.sample_data["face_idx"],
            color="black",
        )

        self._plot_hist()

        # show png
        plt.pause(0.001)
        plt.draw()

    def save_plot(self, output_path):
        plot_path = os.path.join(output_path, "sbm_distribution.png")
        plt.savefig(plot_path)
        plt.close(self.fig)