import numpy as np
from rfcnt import rfcnt, utils # Rainflow module (ftc2)
import os


def __get_script_path():
    return os.path.dirname(__file__)


def example_1():
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        import seaborn as sns
    except ImportError as err:
        print("This example requires modules 'pandas', 'matplotlib' and 'seaborn!'")
        raise err

    data = pd.read_csv(
        os.path.join(__get_script_path(), "long_series.csv"), header=None)
    data = data.to_numpy().squeeze()

    class_count = 50
    class_range = data.max() - data.min()
    class_width = class_range / (class_count - 1)
    class_offset = data.min() - class_width / 2

    res = rfcnt.rfc(
        data, class_count=class_count,
        class_offset=class_offset,
        class_width=class_width,
        hysteresis=class_width,
        use_HCM=0,
        use_ASTM=0,
        #spread_damage=0,   # RFC_SD_HALF_23
        spread_damage=8,   # RFC_SD_TRANSIENT_23c
        #residual_method=0, # RFC_RES_NONE
        residual_method=7, # RFC_RES_REPEATED
        wl={"sd": 1e3, "nd": 1e7, "k": 5})

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(nrows=3, ncols=2, width_ratios=[1, 2])
    ax1 = fig.add_subplot(gs[0, 0])
    sns.heatmap(res["rfm"], cmap="YlOrRd", ax=ax1)
    ax1.invert_yaxis()
    plt.grid(which="both")
    plt.xlabel("Class # (to)")
    plt.ylabel("Class # (from)")
    r = utils.rpplot_prepare(sa=res["rp"][:, 0] / 2, counts=res["rp"][:, 1])
    ax2 = fig.add_subplot(gs[0, 1])
    sns.lineplot(
        x=r["counts"].cumsum(),
        y=r["sa"],
        drawstyle='steps-pre', ci=None, ax=ax2)
    plt.xscale("log")
    plt.ylim(bottom=0, top=2000)
    plt.xlim(left=0.9)
    plt.grid(which="both")
    plt.xlabel("N (log) [1]")
    plt.ylabel("$S_a$")
    ax3 = fig.add_subplot(gs[1, :])
    sns.lineplot(x=np.arange(len(res["dh"])), y=res["dh"].cumsum(), ax=ax3)
    plt.grid(which="both")
    plt.xlabel("Sample #")
    plt.ylabel("Damage (cumulative)")
    ax4 = fig.add_subplot(gs[2, :])
    sns.lineplot(x=np.arange(len(data)), y=data, ax=ax4)
    plt.grid(which="both")
    plt.xlabel("Sample #")
    plt.ylabel("Value")

    fig.tight_layout()
    plt.show()
