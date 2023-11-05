#!/usr/bin/python3

import os
import re

import matplotlib.pyplot as plt
import numpy as np
from colorama import Fore, Style, init as color_init
from tqdm import tqdm

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

LINE_PAT = re.compile(
    r"(\d+) app: (\d+) cumu_bytes: (\d+) cumu_rbs: (\d+) hol_delay: (\d+) "
    r"user: (\d+) slice: (\d+)"
)


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
        # Sacrifice some time for accurate profiling of tqdm
        lines = [line.strip() for line in fin]
        n_lines = len(lines)

    with tqdm(
        total=n_lines,
        ascii=True,
        desc=f"{Fore.MAGENTA + fname + Style.RESET_ALL:<55}",
    ) as pbar:
        # 1: ts, 2: app, 3: cumu_bytes, 4: cumu_rbs, 5: hol_delay, 6: user, 7: slice
        for i, line in enumerate(lines):
            pbar.update(1)
            # Preprocess to avoid regex matching for every line
            words = line.split()
            if len(words) != 13:
                continue
            # Regex matching
            mat = re.match(LINE_PAT, line.strip())
            if mat is None:
                continue
            if int(mat.group(1)) > end_ts:
                pbar.update(n_lines - i - 1)
                break
            if int(mat.group(1)) > begin_ts:
                flow, sid = int(mat.group(2)), int(mat.group(7))
                cumu_rbs[flow] = int(mat.group(4)) / (end_ts / 1000)
                cumu_bytes[flow] = int(mat.group(3)) / (end_ts / 1000)
                flow_to_slice[flow] = sid

    for fid in range(N_USERS):
        sid = flow_to_slice[fid]
        if sid == -1:
            continue
        slice_cumu_bytes[sid] += cumu_bytes[fid]
        slice_cumu_rbs[sid] += cumu_rbs[fid]

    return slice_cumu_rbs, slice_cumu_bytes


def get_throughput_perslice(dname):
    """Returns RB/s, throughput, and sum throughput."""
    ratio = 8 / (1000 * 1000)
    avg_rbs = {k: [0] * N_SLICES for k in BACKLOGGED_MAPPING}
    avg_throughput = {k: [0] * N_SLICES for k in BACKLOGGED_MAPPING}
    sum_throughput = {k: [] for k in BACKLOGGED_MAPPING}

    for i in range(0, TIMES):
        for k, v in BACKLOGGED_MAPPING.items():
            rbs, throughput = get_cumubytes(f"{dname}/backlogged_{v}_{i}.log")
            for j in range(N_SLICES):
                avg_throughput[k][j] += throughput[j] * ratio
                avg_rbs[k][j] += rbs[j]
            sum_throughput[k].append(sum(throughput) * ratio)

    return avg_rbs, avg_throughput, sum_throughput


def plot_together():
    avg_rbs, avg_throughput, sum_throughput = get_throughput_perslice(INPUT_DIR)
    fig, ax = plt.subplots(
        ncols=3, figsize=(26, 7), gridspec_kw={"width_ratios": [2, 2, 1]}
    )

    # Plot 1 & 2: per-slice throughput and per-slice RB/s
    slice_x_array = np.arange(1, N_SLICES + 1)
    plot_config = {
        "bw": {
            "ydata": avg_throughput,
            "ylabel": "Throughput (Mbps)",
            "ylim": {"bottom": 30, "top": 75},
            "title": "(a) Per-slice throughput",
        },
        "rbs": {
            "ydata": avg_rbs,
            "ylabel": "Resource blocks per second (/s)",
            "ylim": {"bottom": 75500, "top": 76750},
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
    bar_x_array = np.arange(0, len(BACKLOGGED_MAPPING), 1)
    bw_array = [np.mean(sum_throughput[k]) for k in BACKLOGGED_MAPPING]
    bwerr_array = [np.std(sum_throughput[k]) for k in BACKLOGGED_MAPPING]

    ax[2].grid(axis="y", alpha=0.4)
    barlist = ax[2].bar(bar_x_array, bw_array, width=0.3)
    for bar, color in zip(barlist, COLORS):
        bar.set_color(color)
    ax[2].errorbar(
        bar_x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black"
    )
    ax[2].set_xlim(bar_x_array[0] - 0.5, bar_x_array[-1] + 0.5)
    ax[2].set_ylim(0, 600)

    ax[2].set_xticks(bar_x_array)
    ax[2].set_xticklabels(BACKLOGGED_MAPPING.keys())
    ax[2].tick_params(axis="x", labelsize=BASE_FONTSIZE + 2)
    ax[2].tick_params(axis="y", labelsize=BASE_FONTSIZE)
    ax[2].set_ylabel("Throughput (Mbps)", fontsize=BASE_FONTSIZE + 2)
    ax[2].set_title("(c) Sum total throughput", y=-0.2, fontsize=BASE_FONTSIZE + 2)

    plt.tight_layout()
    location = f"sameweight.{OUTPUT_TYPE}"
    fig.savefig(location)
    return os.path.abspath(location)


if __name__ == "__main__":
    print()
    color_init(autoreset=True)
    location = plot_together()

    print()
    print(Fore.GREEN + Style.BRIGHT + "Plot available at:", end=" ")
    print(location)
    print()
