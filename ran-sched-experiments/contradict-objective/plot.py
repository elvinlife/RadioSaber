#!/usr/bin/python3

import argparse
import os
import re
import threading
import time
import traceback
from datetime import timedelta
from itertools import cycle, product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from colorama import Fore, Style, init as color_init

INPUT_DIR = "logs"
OUTPUT_DIR = "results"
OUTPUT_TYPE = "png"
TIMES, N_USERS, N_SLICES = 3, 50, 10
BEGIN_TS, END_TS = 0, 10000

VSPACE = 30  # Verbose space, at least 15

AVAIL_FLOW_TYPES = ["vid", "if", "bf"]
AVAIL_INTRA_ALGS = ["mt", "pf", "mlwdf"]

LINE_STYLES = ["X--", "D--", "o--"]
COLORS = ["brown", "dimgrey", "cornflowerblue"]
BASE_FONTSIZE = 18

LOG_MAPPING = {
    "MT": "mt",
    "PF": "pf",
    "M-LWDF": "mlwdf",
}

LINE_PAT = re.compile(
    r"(\d+) app: (\d+) cumu_bytes: (\d+) cumu_rbs: (\d+) hol_delay: (.+) "
    r"user: (\d+) slice: (\d+)"
)

LINE_KEYS = ["ts", "app", "cumu_bytes", "cumu_rbs", "hol_delay", "user", "slice"]

IPFLOW_START_PAT = re.compile(r"ipflow start app: (\d+) flow: (\d+) flowsize: (\d+)")
IPFLOW_END_PAT = re.compile(
    r"ipflow end app: (\d+) flow: (\d+) fct: (.+) flowsize: (\d+) priority: (\d+)"
)

IPFLOW_KEYS = ["fct", "flowsize", "priority"]
IPFLOW_KEYS_SHARED = ["app", "flow"]


def spinner(*, event, prefix, suffix):
    """A spinner in the terminal.

    Parameters
    ----------
    event : threading.Event
        The event used to stop the spinner, done via `event.set()`.
    prefix : str
        The prefix (no space between prefix and spinner).
    suffix : str
        The suffix (no space between spinner and suffix).
    """
    for frame in cycle(r"-\|/"):
        print(f"\r{prefix}{frame}{suffix}", end="", flush=True)
        time.sleep(0.2)
        if event.is_set():
            break


def adjust_ylim(ax, lowerbound=0, upperbound=np.inf):
    """Expand one times the automatic range around the center if possible."""
    cur_ymin, cur_ymax = ax.get_ylim()
    y_adjust = (cur_ymax - cur_ymin) / 2
    ax.set_ylim(
        max(lowerbound, cur_ymin - y_adjust), min(upperbound, cur_ymax + y_adjust)
    )


def get_cdf(data, ratio=0, upperbound=1e5):
    data.sort()
    length = len(data)
    x_array, y_array = [], []
    for i, item in enumerate(data):
        if i / length >= ratio:
            if item > upperbound:
                break
            x_array.append(item)
            y_array.append(i / length)
    return x_array, y_array


def process_log(fname):
    """Process log data.

    Returns
    -------
    slice_cumu_rbs : list of length N_SLICES
        Per-slice cumulative resource blocks per second.
    slice_cumu_bytes : list of length N_SLICES
        Per-slice throughput (bytes per second).
    hol_array : list
        List of head-of-line delays.
    fct_array : list
        List of flow completion times.
    """
    fname_display = Fore.MAGENTA + fname + Style.RESET_ALL

    # Create a spinner for the reading process
    start = time.time()
    event = threading.Event()
    spinner_th = threading.Thread(
        target=spinner,
        kwargs={
            "event": event,
            "prefix": "Reading... ",
            "suffix": " " * (VSPACE - 12) + fname_display,
        },
    )
    spinner_th.start()

    # Read data from the log
    failed, out_of_ts_range = False, False
    try:
        with open(fname, "r", encoding="utf-8") as fin:
            data, data_ipflow = [], {}
            for rawline in fin:
                words = rawline.split()
                line = rawline.strip()

                # For instance:
                # 12062 app: 49 cumu_bytes: 2083456 cumu_rbs: 30224 hol_delay: 0.002 \
                # user: 49 slice: 9
                if len(words) == 13 and words[0].isdigit():
                    mat = re.match(LINE_PAT, line)
                    if mat is None:
                        raise ValueError(f"regex mismatch (len 13, 0th digit)\n{line}")
                    timestamp = int(mat.group(1))
                    if timestamp >= BEGIN_TS and timestamp <= END_TS:
                        data.append(
                            [
                                timestamp,  # ts
                                int(mat.group(2)),  # app
                                int(mat.group(3)),  # cumu_bytes
                                int(mat.group(4)),  # cumu_rbs
                                float(mat.group(5)),  # hol_delay
                                int(mat.group(6)),  # user
                                int(mat.group(7)),  # slice
                            ]
                        )
                    else:
                        out_of_ts_range = True

                # For instance:
                # ipflow start app: 49 flow: 5 flowsize: 10220
                # ipflow end app: 49 flow: 5 fct: 0.01 flowsize: 10220 priority: 0
                elif words[0] == "ipflow":
                    if words[1] == "start" and not out_of_ts_range:
                        mat = re.match(IPFLOW_START_PAT, line)
                        if mat is None:
                            raise ValueError(f"regex mismatch (ipflow start)\n{line}")
                        key = (
                            int(mat.group(1)),  # app
                            int(mat.group(2)),  # flow
                        )
                        data_ipflow[key] = int(mat.group(3))  # flowsize
                    elif words[1] == "end":
                        mat = re.match(IPFLOW_END_PAT, line)
                        if mat is None:
                            raise ValueError(f"regex mismatch (ipflow end)\n{line}")
                        key = (
                            int(mat.group(1)),  # app
                            int(mat.group(2)),  # flow
                        )
                        if key not in data_ipflow:
                            if not out_of_ts_range:
                                raise ValueError(f"ipflow end without start\n{line}")
                            continue  # No start because out of ts range
                        flowsize = int(mat.group(4))  # flowsize
                        if data_ipflow[key] != flowsize:
                            raise ValueError(f"ipflow flowsize mismatch\n{line}")
                        data_ipflow[key] = [
                            float(mat.group(3)),  # fct
                            flowsize,  # flowsize
                            int(mat.group(5)),  # priority
                        ]

        # Main DataFrame
        df = pd.DataFrame(data, columns=LINE_KEYS)

        # IP flow DataFrame (ignoring non-ended IP flows)
        data_ipflow = {k: v for k, v in data_ipflow.items() if not isinstance(v, int)}
        if data_ipflow:
            df_ipflow = pd.DataFrame(data_ipflow).T
            df_ipflow.columns = IPFLOW_KEYS
            df_ipflow = df_ipflow.reset_index(names=IPFLOW_KEYS_SHARED)
        else:
            df_ipflow = pd.DataFrame(columns=IPFLOW_KEYS_SHARED + IPFLOW_KEYS)
    except:  # This is in order to correctly stop the spinner
        failed = True
        print("\r", end="", flush=True)
        traceback.print_exc()

    # Set the event to stop the spinner
    event.set()
    if failed:
        exit(1)
    print("\r", end="", flush=True)

    # Head-of-line delays
    hol_array = df["hol_delay"].tolist()

    # Flow completion times
    fct_array = df_ipflow["fct"].tolist()

    # Per-slice throughput and RB/s
    # df has cumulative data, which needs to be grouped and taken maximum of
    slice_cumu_rbs, slice_cumu_bytes = [], []
    grp = df.groupby(by=["app", "slice"])[["cumu_rbs", "cumu_bytes"]]
    df_cumu = grp.max().reset_index()
    if df_cumu.duplicated(subset="app").any():
        raise ValueError("impossible: same app in different slices")
    for sid in range(N_SLICES):
        df_slice = df_cumu[df_cumu["slice"] == sid]
        slice_cumu_rbs.append(df_slice["cumu_rbs"].sum() / (END_TS / 1000))
        slice_cumu_bytes.append(df_slice["cumu_bytes"].sum() / (END_TS / 1000))

    # Verbosity
    td = timedelta(seconds=int(time.time() - start))
    time_display = f"Done in {td}"
    print(
        Fore.GREEN
        + time_display
        + Style.RESET_ALL
        + " " * (VSPACE - len(time_display))
        + fname_display
    )

    return slice_cumu_rbs, slice_cumu_bytes, hol_array, fct_array


def plot_together(enterprise_alg, flow_type):
    ratio = 8 / (1000 * 1000)  # bytes per second to Mbps
    avg_rbs = {k: [0] * N_SLICES for k in LOG_MAPPING}
    avg_throughput = {k: [0] * N_SLICES for k in LOG_MAPPING}
    sum_throughput = {k: [] for k in LOG_MAPPING}
    all_hols = {k: [] for k in LOG_MAPPING}
    all_fcts = {k: [] for k in LOG_MAPPING}

    for i in range(0, TIMES):
        for k, v in LOG_MAPPING.items():
            slice_rbs, slice_throughput, hols, fcts = process_log(
                f"{INPUT_DIR}/{enterprise_alg}-{flow_type}-{v}-{i}.log"
            )
            for j in range(N_SLICES):
                avg_throughput[k][j] += slice_throughput[j] * ratio
                avg_rbs[k][j] += slice_rbs[j]
            sum_throughput[k].append(sum(slice_throughput) * ratio)
            all_hols[k].append(hols)
            all_fcts[k].append(fcts)

    ### FIRST OUTPUT IMAGE ###

    fig, ax = plt.subplots(
        ncols=3, figsize=(25, 7), gridspec_kw={"width_ratios": [2, 2, 1]}
    )

    # Plot 1 & 2: per-slice throughput and per-slice RB/s
    slice_x_array = np.arange(1, N_SLICES + 1)
    plot_config = {
        "bw": {
            "ydata": avg_throughput,
            "ylabel": "Throughput (Mbps)",
            "title": "(a) Per-slice throughput",
        },
        "rbs": {
            "ydata": avg_rbs,
            "ylabel": "Resource blocks per second (/s)",
            "title": "(b) Per-slice RBs per second",
        },
    }

    for i, config in enumerate(plot_config.values()):
        for k, sty, color in zip(LOG_MAPPING, LINE_STYLES, COLORS):
            ax[i].plot(
                slice_x_array,
                config["ydata"][k],
                sty,
                color=color,
                markersize=8,
                label=k,
            )
        adjust_ylim(ax[i])

        ax[i].set_xlabel("Slice ID", fontsize=BASE_FONTSIZE + 2)
        ax[i].set_ylabel(config["ylabel"], fontsize=BASE_FONTSIZE + 2)
        ax[i].set_xticks(list(range(1, N_SLICES + 1, 2)))
        ax[i].tick_params(axis="both", labelsize=BASE_FONTSIZE)
        ax[i].grid(axis="y", alpha=0.4)
        ax[i].legend(
            fontsize=BASE_FONTSIZE + 2, frameon=False, ncol=3, loc="upper center"
        )
        ax[i].set_title(config["title"], y=-0.2, fontsize=BASE_FONTSIZE + 2)

    # Plot 3: per-slice throughput and per-slice RB/s
    bar_x_array = np.arange(0, len(LOG_MAPPING), 1)
    bw_array = [np.mean(sum_throughput[k]) for k in LOG_MAPPING]
    bwerr_array = [np.std(sum_throughput[k]) for k in LOG_MAPPING]

    ax[2].grid(axis="y", alpha=0.4)
    barlist = ax[2].bar(bar_x_array, bw_array, width=0.3)
    for bar, color in zip(barlist, COLORS):
        bar.set_color(color)
    ax[2].errorbar(
        bar_x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black"
    )
    ax[2].set_xlim(bar_x_array[0] - 0.5, bar_x_array[-1] + 0.5)

    ax[2].set_xticks(bar_x_array)
    ax[2].set_xticklabels(LOG_MAPPING.keys())
    ax[2].tick_params(axis="x", labelsize=BASE_FONTSIZE + 2)
    ax[2].tick_params(axis="y", labelsize=BASE_FONTSIZE)
    ax[2].set_ylabel("Throughput (Mbps)", fontsize=BASE_FONTSIZE + 2)
    ax[2].set_title("(c) Sum total throughput", y=-0.2, fontsize=BASE_FONTSIZE + 2)

    plt.tight_layout()
    location_tprbs = f"{OUTPUT_DIR}/{enterprise_alg}-{flow_type}-tprbs.{OUTPUT_TYPE}"
    fig.savefig(location_tprbs)

    ### SECOND OUTPUT IMAGE ###

    fig, ax = plt.subplots(
        ncols=2, figsize=(20, 7), gridspec_kw={"width_ratios": [1, 1]}
    )

    plot_config = {
        "hol": {
            "cdfdata": all_hols,
            "xlabel": "Queuing delay (s)",
        },
        "fct": {
            "cdfdata": all_fcts,
            "xlabel": "Flow completion time (s)",
        },
    }

    for i, config in enumerate(plot_config.values()):
        for (k, v), color in zip(config["cdfdata"].items(), COLORS):
            ax[i].plot(
                *get_cdf(v[0], ratio=0, upperbound=5),
                "--",
                label=k,
                color=color,
                linewidth=2,
            )
        ax[i].set_xlabel(config["xlabel"], fontsize=BASE_FONTSIZE + 2)
        ax[i].set_ylabel("Ratio", fontsize=BASE_FONTSIZE + 2)
        ax[i].tick_params(axis="both", labelsize=BASE_FONTSIZE)
        ax[i].grid(axis="both", alpha=0.2)
        ax[i].legend(fontsize=BASE_FONTSIZE + 2, frameon=False, loc="lower right")

    plt.tight_layout()
    location_fcthol = f"{OUTPUT_DIR}/{enterprise_alg}-{flow_type}-fcthol.{OUTPUT_TYPE}"
    fig.savefig(location_fcthol)

    return os.path.abspath(location_tprbs), os.path.abspath(location_fcthol)


if __name__ == "__main__":
    color_init(autoreset=True)

    parser = argparse.ArgumentParser(description="Plotting")
    parser.add_argument(
        "--intra-alg",
        choices=AVAIL_INTRA_ALGS + ["all"],
        nargs="+",
    )
    parser.add_argument(
        "--flow-type",
        choices=AVAIL_FLOW_TYPES + ["all"],
        nargs="+",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="run the test plot (ignores other arguments)",
    )
    args = parser.parse_args()

    if args.test:
        # The configuration for the testing is different
        TIMES, N_USERS, N_SLICES = 3, 204, 20
        suites = [("test", "test")]
    else:
        intra_algs, flow_types = args.intra_alg, args.flow_type
        if intra_algs is None or flow_types is None:
            raise ValueError(
                "--intra-alg and --flow-type must be specified, unless using --test "
                "to run the test suite"
            )
        intra_algs = AVAIL_INTRA_ALGS if "all" in intra_algs else set(intra_algs)
        flow_types = AVAIL_FLOW_TYPES if "all" in flow_types else set(flow_types)
        suites = list(product(intra_algs, flow_types))

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    n_suites = len(suites)
    for i, (enterprise_alg, flow_type) in enumerate(suites):
        print()
        print(
            Fore.BLUE
            + Style.BRIGHT
            + f"[{i + 1}/{n_suites}] Enterprise: {enterprise_alg}; Flow: {flow_type}"
        )
        loc_tprbs, loc_fcthol = plot_together(enterprise_alg, flow_type)
        print()
        for location in [loc_tprbs, loc_fcthol]:
            print(Fore.GREEN + Style.BRIGHT + f"Plot available at: {location}")
        print()
