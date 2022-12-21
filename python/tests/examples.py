import os
import numpy as np
import matplotlib.ticker as ticker
from .. import rfc, utils, ResidualMethod, SDMethod


def __get_script_path():
    return os.path.dirname(__file__)


def example_1():
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError as err:
        print("This example requires modules pandas` and `matplotlib`!")
        raise err

    data = pd.read_csv(
        os.path.join(__get_script_path(), "long_series.csv"), header=None)
    data = data.to_numpy().squeeze()

    class_count = 50
    class_range = data.ptp()
    class_width = class_range / (class_count - 1)
    class_offset = data.min() - class_width / 2

    res = rfc(
        data, class_count=class_count,
        class_offset=class_offset,
        class_width=class_width,
        hysteresis=class_width,
        use_HCM=False,
        use_ASTM=False,
        # spread_damage=SDMethod.HALF_23
        spread_damage=SDMethod.TRANSIENT_23c,
        # residual_method=ResidualMethod.NONE,
        residual_method=ResidualMethod.REPEATED,
        wl={"sd": 1e3, "nd": 1e7, "k": 5})

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(nrows=3, ncols=2, width_ratios=[1, 2])

    # RF matrix
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
    plt.grid(which="both")
    plt.xlabel("Class # (to)")
    plt.ylabel("Class # (from)")

    # Range pairs
    r = utils.rpplot_prepare(sa=res["rp"][:, 0] / 2, counts=res["rp"][:, 1])
    ax2 = fig.add_subplot(gs[0, 1])
    # Stairs plot from right to left
    ax2.plot(
        r["counts"].cumsum(),
        r["sa"],
        ds="steps-post")
    plt.xscale("log")
    plt.ylim(bottom=0, top=2500)
    plt.xlim(left=0.9)
    plt.grid(which="both")
    plt.xlabel("N (log) [1]")
    plt.ylabel("$S_a$")

    # Damage history
    ax3 = fig.add_subplot(gs[1, :])
    ax3.plot(np.arange(len(res["dh"])), res["dh"].cumsum())
    plt.grid(which="both")
    plt.xlabel("Sample #")
    plt.ylabel("Damage (cumulative)")

    # Time series
    ax4 = fig.add_subplot(gs[2, :])
    ax4.plot(np.arange(len(data)), data)
    plt.grid(which="both")
    plt.xlabel("Sample #")
    plt.ylabel("Value")

    fig.tight_layout()
    plt.show()
