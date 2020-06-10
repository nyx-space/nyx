# Imports
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.visualization import time_support
from sys import argv
from os import environ


def add_to_plot(hdlr, epochs, dev, cov, idx, name, label):
    color = ["#648FFF", "#DC267F", "#FFB000"][idx % 3]
    max_len = min(epochs.shape[0], dev[idx, :].shape[0])

    epochs = epochs[:max_len]
    y = dev[idx, :][:max_len].flatten()

    hdlr.plot(epochs, y, c=color,
              label="{} deviation".format(name), linewidth=1)
    if cov is not None:
        covs = cov[:, idx][:max_len].flatten()

        hdlr.plot(epochs, covs, linestyle='--', c=color, linewidth=1,
                  label="{} 1-sigma".format(name))
        hdlr.plot(epochs, -covs, linewidth=1,
                  linestyle='--', c=color, label="_nolegend_")

        ylim_x = 1.2 * max(max(y), max(covs))
    else:
        ylim_x = 1.2 * max(abs(y))

    hdlr.set_ylim([-ylim_x, ylim_x])
    hdlr.set_ylabel(label)
    hdlr.grid()
    hdlr.legend()


def main(estimates, truth):
    time_support()

    matplotlib.use('Qt5Agg')

    data = pd.read_csv(estimates).to_numpy()
    truth_data = pd.read_csv(truth).to_numpy()
    print("File loaded.")

    # Extract data
    epochs = (Time(np.array([(epoch.replace("T", " ")[0:-4]) for epoch in data[:, 0]])
                   ) - Time(data[:, 0][0].replace("T", " ")[0:-4])).to_value('hr')
    state_dev = data[:, 1:7].T
    covar_diags = data[:, 7:13]
    est_state = data[:, 13:].T
    truth_state = truth_data[:, 1:7].T
    max_len = min(est_state.shape[1], truth_state.shape[1])
    delta = est_state[:, :max_len] - truth_state[:, :max_len]
    epochs = epochs[:max_len]

    err_avr = ', '.join(['{:1.1e}'.format(sum(compo) / len(compo))
                         for compo in delta[:, :-1]])
    title = "Delta avr. " + err_avr
    print(title)

    node_name = environ.get("PLOT")
    if node_name is None:
        node_name = ""

    node_name += "\n"
    node_name += title

    fig = plt.figure(figsize=(15, 10))
    fig.suptitle(node_name)
    # Pos X states + sigmas
    ax1 = fig.add_subplot(3, 2, 1)
    add_to_plot(ax1, epochs, state_dev, covar_diags,
                0, "x", "Position deviation [km]")

    # Pos Y states + sigmas
    ax2 = fig.add_subplot(3, 2, 3)
    add_to_plot(ax2, epochs, state_dev, covar_diags,
                1, "y", "Position deviation [km]")

    # Pos Z states + sigmas
    ax3 = fig.add_subplot(3, 2, 5)
    add_to_plot(ax3, epochs, state_dev, covar_diags,
                2, "z", "Position deviation [km]")

    # Vel X states + sigmas
    ax4 = fig.add_subplot(3, 2, 2)
    add_to_plot(ax4, epochs, state_dev, covar_diags,
                3, "vx", "Velocity deviation [km/s]")
    # Vel Y states + sigmas
    ax5 = fig.add_subplot(3, 2, 4)
    add_to_plot(ax5, epochs, state_dev, covar_diags,
                4, "vy", "Velocity deviation [km/s]")
    # Vel Z states + sigmas
    ax6 = fig.add_subplot(3, 2, 6)
    add_to_plot(ax6, epochs, state_dev, covar_diags,
                5, "vz", "Velocity deviation [km/s]")

    # Now plot the difference between the estimated state and the real state
    fig = plt.figure(figsize=(15, 10))
    fig.suptitle(node_name)
    # Pos X states + sigmas
    ax1 = fig.add_subplot(3, 2, 1)
    add_to_plot(ax1, epochs, delta, None,
                0, "x", "Position delta [km]")

    # Pos Y states + sigmas
    ax2 = fig.add_subplot(3, 2, 3)
    add_to_plot(ax2, epochs, delta, None,
                1, "y", "Position delta [km]")

    # Pos Z states + sigmas
    ax3 = fig.add_subplot(3, 2, 5)
    add_to_plot(ax3, epochs, delta, None,
                2, "z", "Position delta [km]")

    # Vel X states + sigmas
    ax4 = fig.add_subplot(3, 2, 2)
    add_to_plot(ax4, epochs, delta, None,
                3, "vx", "Velocity delta [km/s]")
    # Vel Y states + sigmas
    ax5 = fig.add_subplot(3, 2, 4)
    add_to_plot(ax5, epochs, delta, None,
                4, "vy", "Velocity delta [km/s]")
    # Vel Z states + sigmas
    ax6 = fig.add_subplot(3, 2, 6)
    add_to_plot(ax6, epochs, delta, None,
                5, "vz", "Velocity delta [km/s]")

    fig.show()
    plt.show()


if __name__ == "__main__":
    if len(argv) == 3:
        est = argv[1]
        trt = argv[2]
    else:
        est = "data/estimates.csv"
        trt = "data/truth.csv"
    main(est, trt)
