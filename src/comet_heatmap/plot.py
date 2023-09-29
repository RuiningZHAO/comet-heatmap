import h5py

# NumPy
import numpy as np
# AstroPy
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table
# Plot
from matplotlib import rc, cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from . import conf

# Set plot parameters
rc('axes', **{'lw': 1})
rc('font', **{'family': 'serif'})
rc('text', usetex=True)


def heatmap(data, row_labels, col_labels, colormap, ax=None, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    cmap = cm.get_cmap(colormap)
    cmap.set_over("w")
    cmap.set_bad("lightgrey")
    im = ax.imshow(data, cmap = cmap, **kwargs)
    ax.grid(which="minor", color="k", linestyle="-", linewidth=1)

    # Show 1/3 x ticks...
    ax.set_xticks(np.arange(data.shape[1])[::3], minor=False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]),      minor=False)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # Let the horizontal axes labeling appear on top.
    ax.tick_params(
        which='major', width=1, top=True, bottom=False, labeltop=True, 
        labelbottom=False)
    ax.tick_params(which="minor", bottom=False, left=False)
    # Set tick labels
    ax.set_xticklabels(
        col_labels[::3], rotation=-30, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(row_labels)

    return im


def plot(save_dir, mag_lim=conf.limiting_magnitude, 
         path_to_library=conf.path_to_library):
    """
    """
    
    # Load
    with h5py.File(os.path.join(path_to_library, "matrices.hdf5", "r") as hf:
        # 1. Catalog
        catalog = Table.read(hf, path="catalog")
        catalog["designation"] = np.char.decode(catalog["designation"].value)
        catalog["reference"]   = np.char.decode(catalog["reference"  ].value)
        # 2. Time & FLI
        FLI = Table.read(hf, path="FLI")
        time_tag = FLI["time_tag"].to_value("iso", subfmt="date") # UTC???
        lunar_illum = FLI["lunar_illum"].value
        # 3. Tmag
        Tmag_mat = hf["Tmag"][:]
        # 4. is high
        is_high_mat = hf["is_high"][:]

    # Sort
    Tmag_mat[~is_high_mat] = np.nan
    idx = np.argsort(np.nanmin(Tmag_mat, axis=0))
    catalog = catalog[idx]
    Tmag_mat = Tmag_mat[:, idx]
    is_high_mat = is_high_mat[:, idx]

    # Filter
    mask_non = catalog["eccentricity"] >= 1
    mask_per = ~mask_non
    
    # Plot
    for mask, colormap, fig_type in zip(
        [mask_non, mask_per], ["winter", "autumn"], ["non", "per"]):

        figsize = (
            0.15 * (Tmag_mat[:, mask].shape[0] + 5.5), 
            0.15 * (Tmag_mat[:, mask].shape[1] + 10)
        )

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        # Heatmap
        im = heatmap(
            Tmag_mat[:, mask].T, catalog["designation"].value[mask], time_tag, 
            colormap=colormap, vmax=mag_lim, ax=ax)
        # Add axes
        divider = make_axes_locatable(ax)
        # Moon phase
        sax = divider.append_axes("bottom", size=0.15, pad=0.05)
        sax.scatter(
            np.arange(time_tag.shape[0]), np.zeros(time_tag.shape[0]) + mask.sum(), 
            s=50, c=lunar_illum, edgecolors="k", linewidths=1, cmap="Greys_r")
        sax.axis("off")
        sax.set_xlim(-0.5, time_tag.shape[0] - 0.5)
        # Colorbar
        cax = divider.append_axes("bottom", size=0.15, pad=0.05)
        cbar = cax.figure.colorbar(
            im, cax=cax, label="Total magnitude", orientation="horizontal")
        cbar.ax.tick_params(which="major", width=1)

        fig.tight_layout()
        plt.savefig(
            os.path.join(save_dir, f"comet_observability_{fig_type}.png"), dpi=144)


def debug():
    """
    """

    plot(mag_lim=20)

    return None


if __name__ == "__main__":
    debug()
