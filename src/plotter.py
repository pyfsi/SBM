from utils import np, sys, os, plt

class Plotter():
    def __init__(self):
        # distribution variables
        self.variables = {"mass": "kg", "position": "m"}#, "position-y": "m", "position-z": "m"}

        # fig and axes
        n_vars = len(self.variables.keys())
        self.fig, self.axs = plt.subplots(n_vars, 2, layout='constrained')
        self.axs[0,0].title.set_text('Successful samples')
        self.axs[0,1].title.set_text('Rejected samples')
        for i, (var, unit) in enumerate(self.variables.items()):
            # ax i,0 settings
            self.axs[i,0].set_ylabel('Count [-]')
            self.axs[i,0].set_xlabel(f'{var} [{unit}]')
            # ax i,1 settings
            self.axs[i,1].set_ylabel('Count [-]')
            self.axs[i,1].set_xlabel(f'{var} [{unit}]')

    def initialize(self, config, inlet_data):
        self.output_path = config["_output_path"]
        # mass bounds
        self.mass_min = float(config["model"]["mass_g"]["min"])
        self.mass_max = float(config["model"]["mass_g"]["max"])

        # position bounds, reduced to 1D using L2 norm
        self.pos_min = 0.0
        self.pos_max = max(np.linalg.norm(inlet_data.min, ord=2), np.linalg.norm(inlet_data.max, ord=2))

        # inlet geometry
        self.inlet_center = inlet_data.center
        self.inlet_normal = inlet_data.normal

        # distribution parameters for mass
        mass_bins, mass_counts = self._calc_bins_and_counts(self.mass_min, self.mass_max, 100)
        pos_bins, pos_counts = self._calc_bins_and_counts(self.pos_min, self.pos_max, 100)

        self.distribution = {}
        self.distribution_rejected = {}
        self.oob_samples = 0 # for bubbles outside of mass bounds

        # set values
        self.distribution["mass"] = (np.copy(mass_bins), np.copy(mass_counts))
        self.distribution_rejected["mass"] = (np.copy(mass_bins), np.copy(mass_counts))
        self.distribution["position"] = (np.copy(pos_bins), np.copy(pos_counts))
        self.distribution_rejected["position"] = (np.copy(pos_bins), np.copy(pos_counts))

    def plot(self, samples, rejected_samples):
        self.update(samples, rejected_samples)
        self.draw()

    def update(self, samples, rejected_samples):
        # ignore and count mass samples out of bounds
        samples_pass = {"mass":[], "position":[]}
        for m, x in zip(samples["mass"], samples["position"]):
            if (m<self.mass_min) | (m>self.mass_max):
                self.oob_samples += 1
            else:
                samples_pass["mass"].append(m)
                proj_point = self._project_points(x)
                samples_pass["position"].append(np.linalg.norm(proj_point))

        # process rejected samples
        samples_fail = {"mass":[], "position":[]}
        for m, x in zip(rejected_samples["mass"], rejected_samples["position"]):
            samples_fail["mass"].append(m)
            proj_point = self._project_points(x)
            samples_fail["position"].append(np.linalg.norm(proj_point))

        # succcessful samples
        pos_idx = np.digitize(samples_pass["position"], self.distribution["position"][0]) - 1
        np.add.at(self.distribution["position"][1], pos_idx, 1)
        mass_idx = np.digitize(samples_pass["mass"], self.distribution["mass"][0]) - 1
        np.add.at(self.distribution["mass"][1], mass_idx, 1)

        # rejected samples
        pos_fail_idx = np.digitize(samples_fail["position"], self.distribution_rejected["position"][0]) - 1
        id_min, id_max = 0, len(self.distribution_rejected["position"][1])
        # pos_fail_idx = pos_fail_idx[(id_min<=pos_fail_idx)&(pos_fail_idx<id_max)]
        np.add.at(self.distribution_rejected["position"][1], pos_fail_idx, 1)
        mass_fail_idx = np.digitize(samples_fail["mass"], self.distribution_rejected["mass"][0]) - 1
        np.add.at(self.distribution_rejected["mass"][1], mass_fail_idx, 1)

    def draw(self):
        # remove previous figures
        self._clear_figures()
        self._plot_histograms()

        # show png
        plt.pause(0.001)
        plt.draw()

    def hold_plot(self):
        plt.show()

    def save_plot(self):
        plot_path = os.path.join(self.output_path, "sbm_distribution.png")
        plt.savefig(plot_path)
        plt.close(self.fig)

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

    def _plot_histograms(self):
        for row, var in enumerate(self.variables.keys()):
            distr = self.distribution[var]
            self._draw_histogram(distr, row, 0, color="green")
            distr_rejected =  self.distribution_rejected[var]
            self._draw_histogram(distr_rejected, row, 1, color="red")

    def _draw_histogram(self, distr, row, col, color="green"):
        # plot settings template
        box_props = dict(boxstyle='round', edgecolor='black', facecolor="white", alpha=0.9)
        bar_kwargs = {"x":None, "height":None,
                        "edgecolor":"black", "facecolor":None, "width":None,}
        txt_kwargs = {"x":0.02, "y":0.93, "s":"",
                        "transform":None, "fontsize":"15","bbox":box_props}

        n_bubbles = np.sum(distr[1])
        # calculate width of bar plot assuming same spacing
        width = (distr[0][1]-distr[0][0])
        bar_x = distr[0][:] + 0.5 * width

        # assign to kwargs
        bar_kwargs["x"] = bar_x
        bar_kwargs["height"] = distr[1][:]
        bar_kwargs["facecolor"] = color
        bar_kwargs["width"] = width
        txt_kwargs["transform"] = self.axs[row,col].transAxes
        txt_kwargs["s"] = f"n={int(n_bubbles)}"

        # plot
        self.axs[row,col].bar(**bar_kwargs)
        self.axs[row,col].text(**txt_kwargs)

    def _calc_bins_and_counts(self, min, max, n_bins):
        bin_width = (max - min) / n_bins
        min_id = np.round(min/bin_width)
        max_id = np.round(max/bin_width)
        bins = np.arange(min_id, max_id+1, 1) * bin_width
        counts = np.zeros(n_bins+1)
        return bins, counts

    def _project_points(self, point):
        '''
        Project a point in 3D space onto the inlet plane.
        Args:
            point: point coordinates
            plane_center: coordinate of plane center
            plane_normal: plane normal vector
        Returns:
            vector from plane_center to the point's projection
        '''
        unit_plane_normal = self.inlet_normal / np.linalg.norm(self.inlet_normal)
        relative_position = point - self.inlet_normal
        return relative_position - np.dot(relative_position, unit_plane_normal) * unit_plane_normal