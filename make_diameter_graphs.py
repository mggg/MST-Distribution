import json
import jsonlines
import glob
import matplotlib.pyplot as plt
from fractions import Fraction
from pathlib import Path


def make_plot(data, save_path):
    diameters = []
    probabilities = []
    frequencies = []

    for diameter, prob_dict in data.items():
        for prob, freq in prob_dict.items():
            diameters.append(diameter)
            probabilities.append(Fraction(prob))
            frequencies.append(freq)

    max_freq = max(frequencies)
    sum_freq = sum(frequencies)

    frequencies = [freq / max_freq * 300 for freq in frequencies]

    plt.scatter(
        diameters,
        probabilities,
        s=frequencies,
        alpha=0.5,
        c="blue",
        edgecolors="w",
        linewidth=0.5,
    )

    plt.xlabel("Diameter")
    plt.ylabel("Probability")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":

    for file_name in glob.glob("./enumerations/*"):
        with jsonlines.open(file_name, "r") as reader:
            plot_data = {}
            # {diameter: {prob: freq}}
            for obj in reader:
                diameter = obj["diameter"]
                prob = obj["probability"]
                plot_data[diameter] = plot_data.get(diameter, {})
                plot_data[diameter][prob] = plot_data[diameter].get(prob, 0) + 1

            file_name = Path.joinpath(
                Path("./diameter_vs_probability_plots"),
                Path(Path(file_name).stem + "_diameter_vs_probability.png"),
            )
            make_plot(plot_data, file_name)
