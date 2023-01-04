import os
import numpy as np
from .. import ResidualMethod, SDMethod, rfc, utils, ArrayLike


def _get_script_path() -> str:
    """
    Get the directory path of the current script.

    Returns
    -------
    str
        The directory path where the current script is located.
    """
    return os.path.dirname(os.path.abspath(__file__))


def example_1():
    """
    Example function to demonstrate rainflow cycle counting and data visualization.

    This function reads a time series data from a CSV file, performs rainflow cycle counting,
    and generates plots for the rainflow matrix, range pairs, damage history, and the original
    time series data.

    Raises
    ------
    ImportError
        If the required modules 'pandas' and 'matplotlib' are not installed.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        import pandas as pd
        from matplotlib.gridspec import GridSpec
    except ImportError as err:
        print("This example requires modules 'pandas' and 'matplotlib'!")
        raise err

    # Load data from CSV
    data = pd.read_csv(os.path.join(_get_script_path(), "long_series.csv"), header=None)
    data = data.to_numpy().squeeze()

    # Set parameters for rainflow counting
    class_count = 50
    class_range = data.ptp()
    class_width = class_range / (class_count - 1)
    class_offset = data.min() - class_width / 2

    # Perform rainflow counting
    res = rfc(
        data,
        class_count=class_count,
        class_offset=class_offset,
        class_width=class_width,
        hysteresis=class_width,
        use_HCM=False,
        use_ASTM=False,
        spread_damage=SDMethod.TRANSIENT_23c,
        residual_method=ResidualMethod.REPEATED,
        wl={"sd": 1e3, "nd": 1e7, "k": 5},
    )

    # Create a figure and gridspec layout
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(nrows=3, ncols=2, width_ratios=[1, 2])

    # Rainflow matrix plot
    ax1 = fig.add_subplot(gs[0, 0])
    im = ax1.imshow(res["rfm"], cmap="YlOrRd", aspect=0.7)
    cb = plt.colorbar(im, ax=ax1, label="Counts")
    cb.ax.tick_params(labelsize="x-small")
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax1.xaxis.set_ticks_position("top")
    ax1.xaxis.set_label_position("top")
    ax1.tick_params(axis="x", labelsize="x-small", labelrotation=90)
    ax1.tick_params(axis="y", labelsize="x-small")
    ax1.grid(which="both")
    ax1.set_xlabel("Class # (to)")
    ax1.set_ylabel("Class # (from)")

    # Range pairs plot
    r = utils.rpplot_prepare(sa=res["rp"][:, 0] / 2, counts=res["rp"][:, 1])
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(r["counts"].cumsum(), r["sa"], drawstyle="steps-post")
    ax2.set_xscale("log")
    ax2.set_ylim(bottom=0, top=2500)
    ax2.set_xlim(left=0.9)
    ax2.grid(which="both")
    ax2.set_xlabel("N (log) [1]")
    ax2.set_ylabel("$S_a$")

    # Damage history plot
    ax3 = fig.add_subplot(gs[1, :])
    ax3.plot(np.arange(len(res["dh"])), res["dh"].cumsum())
    ax3.grid(which="both")
    ax3.set_xlabel("Sample #")
    ax3.set_ylabel("Damage (cumulative)")

    # Time series plot
    ax4 = fig.add_subplot(gs[2, :])
    ax4.plot(np.arange(len(data)), data)
    ax4.grid(which="both")
    ax4.set_xlabel("Sample #")
    ax4.set_ylabel("Value")

    fig.tight_layout()
    plt.show()
