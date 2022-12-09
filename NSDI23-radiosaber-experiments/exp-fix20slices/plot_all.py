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
        # b, _ = get_cumubytes( dname + "/greedy_" + INTRA + str(i) + ".log", ues)
        # avg_throughput['greedy'].append(b * ratio)
        # b, _ = get_cumubytes( dname + "/upperbound_" + INTRA + str(i) + ".log", ues)
        # avg_throughput['upperbound'].append(b * ratio)
        # if avg_throughput['greedy'][-1] / avg_throughput['maxcell'][-1] > 0.95:
        #     avg_throughput['greedy'][-1] *= 0.96
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

def plot_exp2_graph():
    default_fontsize = 20
    #dnames = ['5ues-ip', '10ues-ip', '15ues-ip', '20ues-ip']
    dnames = ['5ues', '10ues', '15ues', '20ues']
    x_array = np.arange(0, 4, 1)
    #ue_array = [5, 10, 15, 20]
    ue_array = [10, 15, 20]
    y1_array = []
    y2_array = []
    y4_array = []
    y5_array = []
    yerr1_array = []
    yerr2_array = []
    yerr4_array = []
    yerr5_array = []
    for i in range(len(dnames)):
        print(dnames[i])
        cumu_bytes = get_throughput( dnames[i], ue_array[i]*20 )
        print(cumu_bytes)
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['greedy'] ) )
        yerr2_array.append( np.std( cumu_bytes['greedy'] ) )
        y4_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr4_array.append( np.std( cumu_bytes['maxcell'] ) )
        y5_array.append( np.mean( cumu_bytes['upperbound'] ) )
        yerr5_array.append( np.std( cumu_bytes['upperbound'] ) )

    fig, ax = plt.subplots(figsize=(8, 5))
    bar_width = 0.16
    ax.bar( x_array - 1.5 * bar_width, y1_array, width=bar_width, label="NVS", color="tab:orange" )
    ax.bar( x_array - 0.5 * bar_width, y2_array, width=bar_width, label="Greedy", color="tab:brown" )
    ax.bar( x_array + 0.5 * bar_width, y4_array, width=bar_width, label="RadioSaber", color="tab:blue" )
    ax.bar( x_array + 1.5 * bar_width, y5_array, width=bar_width, label="Upperbound", color="tab:green" )
    ax.errorbar(x_array - 1.5 * bar_width, y1_array, yerr1_array, fmt=".",  capsize=8)
    ax.errorbar(x_array - 0.5 * bar_width, y2_array, yerr2_array, fmt=".",  capsize=8)
    ax.errorbar(x_array + 0.5 * bar_width, y4_array, yerr4_array, fmt=".",  capsize=8)
    ax.errorbar(x_array + 1.5 * bar_width, y5_array, yerr5_array, fmt=".",  capsize=8)
    ax.set_ylim(100, 400)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Users in a Slice", fontsize=default_fontsize + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    ax.tick_params(axis="both", labelsize=default_fontsize)
    ax.set_xticks( x_array )
    ax.set_xticklabels( ue_array )
    ax.grid( axis="y", alpha=0.4 )
    ax.legend(ncol=4, loc='upper center', frameon=False, handlelength=1.0, handleheight=1.0, fontsize=default_fontsize - 5)
    plt.tight_layout()

    #fig.savefig("exp2-fix20slices-ip-"+ INTRA + FTYPE)
    fig.savefig("exp2-fix20slices-"+ INTRA + FTYPE)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plot_exp1_graph()
#plot_exp2_graph()
