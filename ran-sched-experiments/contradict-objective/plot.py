#!/usr/bin/python3


import matplotlib.pyplot as plt
import numpy as np

INPUT_DIR = "results"
OUTPUT_TYPE = "png"
TIMES, N_USERS, N_SLICES = 3, 204, 20

LINE_STYLES = ["X--", "D--", "o--"]
COLORS = ["brown", "dimgrey", "cornflowerblue"]
BASE_FONTSIZE = 18

BACKLOGGED_MAPPING = {
    "MT": "max_throughput",
    "PF": "proportional_fairness",
    "M-LWDF": "m_lwdf",
}


def get_cdf(data, ratio=0):
    data.sort()
    x_array = []
    y_array = []
    for i, item in enumerate(data):
        if i / len(data) >= ratio:
            x_array.append(item)
            y_array.append(i / len(data))
    return x_array, y_array


def get_cumubytes(fname):
    """Get the per-second cumulative sent bytes."""
    begin_ts = 0
    end_ts = 10000

    cumu_bytes = [0] * N_USERS
    cumu_rbs = [0] * N_USERS
    flow_to_slice = [-1] * N_USERS
    slice_cumu_bytes = [0] * N_SLICES
    slice_cumu_rbs = [0] * N_SLICES

    with open(fname, "r", encoding="utf-8") as fin:
        for line in fin:
            words = line.split()
            # Check if there are enough elements in the 'words' list to prevent 'IndexError'
            if len(words) < 13:
                # print(len(words))
                continue  # Skip lines that don't have enough elements
            if not words[0].isdigit():
                continue
            if int(words[0]) > end_ts:
                break
            if int(words[0]) > begin_ts:
                # Parse the log entry correctly based on the labels in the format
                flow = int(words[words.index("app:") + 1])
                sid = int(words[words.index("slice:") + 1])
                cumu_rbs[flow] = int(words[words.index("cumu_rbs:") + 1]) / (
                    end_ts / 1000
                )
                cumu_bytes[flow] = int(words[words.index("cumu_bytes:") + 1]) / (
                    end_ts / 1000
                )
                flow_to_slice[flow] = sid

    for fid in range(N_USERS):
        sid = flow_to_slice[fid]
        if sid == -1:
            continue
        slice_cumu_bytes[sid] += cumu_bytes[fid]
        slice_cumu_rbs[sid] += cumu_rbs[fid]

    return slice_cumu_rbs, slice_cumu_bytes


def get_throughput_perslice(dname):
    ratio = 8 / (1000 * 1000)
    avg_throughput = {k: [0] * N_SLICES for k in BACKLOGGED_MAPPING}
    avg_rbs = {k: [0] * N_SLICES for k in BACKLOGGED_MAPPING}

    for i in range(0, TIMES):
        for k, v in BACKLOGGED_MAPPING.items():
            rbs, throughput = get_cumubytes(f"{dname}/backlogged_{v}_{i}.log")
            for j in range(N_SLICES):
                avg_throughput[k][j] += throughput[j] * ratio
                avg_rbs[k][j] += rbs[j]

    return avg_rbs, avg_throughput


def get_throughput_total(dname):
    ratio = 8 / (1000 * 1000)
    avg_throughput = {k: [] for k in BACKLOGGED_MAPPING}

    for i in range(0, TIMES):
        for k, v in BACKLOGGED_MAPPING.items():
            _, throguhput = get_cumubytes(f"{dname}/backlogged_{v}_{i}.log")
            avg_throughput[k].append(sum(throguhput) * ratio)

    return avg_throughput


def plot_together():
    fig, ax = plt.subplots(
        ncols=3, figsize=(26, 7), gridspec_kw={"width_ratios": [2, 2, 1]}
    )

    # Plot 1 & 2: per-slice throughput and per-slice RB/s
    avg_rbs, avg_throughput = get_throughput_perslice(INPUT_DIR)
    slice_x_array = np.arange(1, N_SLICES + 1)

    plot_config = {
        "bw": {
            "ydata": avg_throughput,
            "ylabel": "Throughput (Mbps)",
            "ylim": {"bottom": 0, "top": 100},
            "title": "(a) Per-slice throughput",
        },
        "rbs": {
            "ydata": avg_rbs,
            "ylabel": "Resource blocks per second (/s)",
            "ylim": {"bottom": 75000, "top": 80000},
            "title": "(b) Per-slice RBs per second",
        },
    }

    for i, (k, config) in enumerate(plot_config.items()):
        for k, sty, color in zip(BACKLOGGED_MAPPING, LINE_STYLES, COLORS):
            ax[i].plot(
                slice_x_array,
                config["ydata"][k],
                sty,
                color=color,
                markersize=8,
                label=k,
            )

        ax[i].set_xlabel("Slice ID", fontsize=BASE_FONTSIZE + 2)
        ax[i].set_ylabel(config["ylabel"], fontsize=BASE_FONTSIZE + 2)
        ax[i].set_ylim(**config["ylim"])
        ax[i].set_xticks(list(range(1, N_SLICES + 1, 2)))
        ax[i].tick_params(axis="both", labelsize=BASE_FONTSIZE)
        ax[i].grid(axis="y", alpha=0.4)
        ax[i].legend(
            fontsize=BASE_FONTSIZE + 2, frameon=False, ncol=3, loc="upper center"
        )
        ax[i].set_title(config["title"], y=-0.2, fontsize=BASE_FONTSIZE + 2)

    # Plot 3: per-slice throughput and per-slice RB/s
    sum_throughput = get_throughput_total(INPUT_DIR)
    bar_x_array = np.arange(0, len(BACKLOGGED_MAPPING), 1)

    bw_array = [np.mean(sum_throughput[k]) for k in BACKLOGGED_MAPPING]
    bwerr_array = [np.std(sum_throughput[k]) for k in BACKLOGGED_MAPPING]

    ax[2].grid(axis="y", alpha=0.4)
    barlist = ax[2].bar(bar_x_array, bw_array, width=0.2)
    for bar, color in zip(barlist, COLORS):
        bar.set_color(color)
    ax[2].errorbar(
        bar_x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black"
    )
    ax[2].set_xlim(bar_x_array[0] - 0.5, bar_x_array[-1] + 0.5)
    ax[2].set_ylim(0, 700)

    ax[2].set_xticks(bar_x_array)
    ax[2].set_xticklabels(BACKLOGGED_MAPPING.keys())
    ax[2].tick_params(axis="x", labelsize=BASE_FONTSIZE + 2)
    ax[2].tick_params(axis="y", labelsize=BASE_FONTSIZE)
    ax[2].set_ylabel("Throughput (Mbps)", fontsize=BASE_FONTSIZE + 2)
    ax[2].set_title("(c) Sum total throughput", y=-0.2, fontsize=BASE_FONTSIZE + 2)

    plt.tight_layout()
    fig.savefig(f"sameweight.{OUTPUT_TYPE}")


if __name__ == "__main__":
    plot_together()
