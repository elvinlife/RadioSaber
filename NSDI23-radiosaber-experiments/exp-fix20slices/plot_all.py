#!/usr/bin/python3
TIMES=1
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
import numpy as np

INTRA = "pf"
FTYPE=".pdf"
COLORS=["dimgrey", "", "cornflowerblue", ""]

# get the per-second cumulative bytes
def get_cumubytes(fname, n_users):
    end_ts = 10000
    cumu_bytes = [0 for i in range(n_users)]
    cumu_rbs = [0 for i in range(n_users)]
    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if not words[0].isdigit():
                continue
            if int(words[0]) > end_ts:
                break
            flow = int(words[2])
            cumu_rbs[flow] = int( words[6] ) / (end_ts / 1000)
            cumu_bytes[flow] = int( words[4] ) / (end_ts / 1000)

    return sum(cumu_bytes), sum(cumu_rbs)

def get_throughput(dname, ues):
    ratio = 8 / (1000 * 1000)
    avg_throughput = defaultdict(list)
    rbs_dict = defaultdict(list)
    for i in range(TIMES):
        cumu_bytes, cumu_rbs = get_cumubytes( dname + "/nvs_" + INTRA + str(i) + ".log", ues)
        avg_throughput['nvs'].append(cumu_bytes * ratio)
        rbs_dict['nvs'].append(cumu_rbs)
        cumu_bytes, cumu_rbs = get_cumubytes( dname + "/maxcell_" + INTRA + str(i) + ".log", ues)
        avg_throughput['maxcell'].append(cumu_bytes * ratio)
        rbs_dict['maxcell'].append(cumu_rbs)
    return avg_throughput

def plot_exp1_graph():
    default_fontsize = 24
    dnames = ['5ues-ip', '10ues-ip', '15ues-ip', '20ues-ip']
    # dnames = ['5ues', '10ues', '15ues', '20ues']
    x_array = np.arange(0, 4, 1)
    ue_array = [5, 10, 15, 20]
    y1_array = []
    y2_array = []
    yerr1_array = []
    yerr2_array = []
    for i in range(len(dnames)):
        print(dnames[i])
        cumu_bytes = get_throughput( dnames[i], ue_array[i]*20 )
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr2_array.append( np.std( cumu_bytes['maxcell'] ) )

    fig, ax = plt.subplots(figsize=(8, 5))
    bar_width = 0.2
    ax.bar( x_array - 0.5 * bar_width, y1_array, width=bar_width, label="NVS", color=COLORS[0] )
    ax.bar( x_array + 0.5 * bar_width, y2_array, width=bar_width, label="RadioSaber", color=COLORS[2] )
    ax.errorbar(x_array - 0.5 * bar_width, y1_array, yerr1_array, fmt=".", capsize=10, color="black")
    ax.errorbar(x_array + 0.5 * bar_width, y2_array, yerr2_array, fmt=".", capsize=10, color="black")
    ax.set_ylim(0, 300)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Num of UEs per Slice", fontsize=default_fontsize + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    ax.tick_params(axis="both", labelsize=default_fontsize-2)
    ax.set_xticks( x_array )
    ax.set_xticklabels( ue_array )
    ax.grid( axis="y", alpha=0.4 )
    ax.legend(ncol=2, loc='upper center', frameon=False,\
            handlelength=1.0, handleheight=1.0,\
            bbox_to_anchor=(0.5, 1.05),\
            fontsize=default_fontsize + 2)
    plt.tight_layout()

    fig.savefig("exp1-fix20slices-ip-"+ INTRA + FTYPE)
    # fig.savefig("exp1-fix20slices-"+ INTRA + FTYPE)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plot_exp1_graph()
