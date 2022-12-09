#!/usr/bin/python3
TIMES=1
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
INTRA = "pf"
FTYPE=".pdf"
COLORS=["dimgrey", "brown", "cornflowerblue"]

# get the per-second cumulative bytes
def get_cumubytes(fname, n_users):
    begin_ts = 9900
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
            if int(words[0]) > begin_ts:
                flow = int(words[2])
                cumu_rbs[flow] = int( words[6] ) / (end_ts / 1000)
                cumu_bytes[flow] = int( words[4] ) / (end_ts / 1000)

    return sum(cumu_bytes), sum(cumu_rbs)

def get_throughput(dname, ues):
    ratio = 8 / (1000 * 1000)
    cumu_bytes = defaultdict(list)
    for i in range(1, TIMES+1):
        b, _ = get_cumubytes( dname + "/nvs_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['nvs'].append(b * ratio)
        b, _ = get_cumubytes( dname + "/nvs_optimal_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['nvs_optimal'].append(b * ratio)
        b, _ = get_cumubytes( dname + "/maxcell_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['maxcell'].append(b * ratio)
    return cumu_bytes

def plot_exp_graph():
    default_fontsize = 20
    fig, ax = plt.subplots(figsize=(12, 5))
    bar_width = 0.14

    dnames = ['10uesX20slices', '20uesX20slices', '30uesX20slices']
    x_array = np.arange(0, 3, 1)
    ue_array = [10, 20, 30]
    y1_array = []
    y2_array = []
    y3_array = []
    yerr1_array = []
    yerr2_array = []
    yerr3_array = []
    for i in range(len(dnames)):
        cumu_bytes = get_throughput( dnames[i], ue_array[i]*20 )
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['nvs_optimal'] ) )
        yerr2_array.append( np.std( cumu_bytes['nvs_optimal'] ) )
        y3_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr3_array.append( np.std( cumu_bytes['maxcell'] ) )
        print( "nvs_nongreed->nvs: %f nvs_nongreed->maxcell: %f" % (
            np.mean(cumu_bytes['nvs_optimal']) / np.mean(cumu_bytes['nvs']) - 1,
            1 - np.mean(cumu_bytes['nvs_optimal']) / np.mean(cumu_bytes['maxcell']),
            ))

    ax.bar( x_array - bar_width, y1_array, width=bar_width, label="NVS w/ greedy PF", color=COLORS[0] )
    ax.bar( x_array , y2_array, width=bar_width, label="NVS w/ non-greedy PF", color=COLORS[1] )
    ax.bar( x_array + bar_width, y3_array, width=bar_width, label="RadioSaber", color=COLORS[2] )
    ax.errorbar(x_array - bar_width, y1_array, yerr1_array, fmt=".", capsize=15, color="black")
    ax.errorbar(x_array, y2_array, yerr2_array, fmt=".", capsize=15, color="black")
    ax.errorbar(x_array + bar_width, y3_array, yerr3_array, fmt=".", capsize=15, color="black")
    ax.set_ylim(100, 300)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Num of UEs per Slice", fontsize=default_fontsize + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    ax.tick_params(axis="both", labelsize=default_fontsize)
    ax.set_xticks( x_array )
    ax.set_xticklabels( ue_array )
    ax.grid( axis="y", alpha=0.4 )
    ax.legend(ncol=4, loc='upper center', frameon=False, handlelength=1.0, handleheight=1.0, fontsize=default_fontsize )

    # dnames = ['20uesX10slices', '20uesX20slices', '20uesX30slices']
    # x_array = np.arange(0, 3, 1)
    # slice_array = [10, 20, 30]
    # y1_array = []
    # y2_array = []
    # y3_array = []
    # yerr1_array = []
    # yerr2_array = []
    # yerr3_array = []
    # for i in range(len(dnames)):
    #     cumu_bytes = get_throughput( dnames[i], slice_array[i]*20 )
    #     y1_array.append( np.mean( cumu_bytes['nvs'] ) )
    #     yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
    #     y2_array.append( np.mean( cumu_bytes['nvs_optimal'] ) )
    #     yerr2_array.append( np.std( cumu_bytes['nvs_optimal'] ) )
    #     y3_array.append( np.mean( cumu_bytes['maxcell'] ) )
    #     yerr3_array.append( np.std( cumu_bytes['maxcell'] ) )

    # ax[1].bar( x_array - bar_width, y1_array, width=bar_width, label="NVS", color="tab:orange" )
    # ax[1].bar( x_array , y2_array, width=bar_width, label="NVS(Optimal)", color="tab:brown" )
    # ax[1].bar( x_array + bar_width, y3_array, width=bar_width, label="RadioSaber", color="tab:blue" )
    # ax[1].errorbar(x_array - bar_width, y1_array, yerr1_array, fmt=".",  capsize=8)
    # ax[1].errorbar(x_array, y2_array, yerr2_array, fmt=".",  capsize=8)
    # ax[1].errorbar(x_array + bar_width, y3_array, yerr3_array, fmt=".",  capsize=8)
    # ax[1].set_ylim(100, 300)
    # ax[1].spines["top"].set_visible(False)
    # ax[1].spines["right"].set_visible(False)
    # ax[1].set_xlabel("Users in a Slice", fontsize=default_fontsize + 4)
    # ax[1].set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    # ax[1].tick_params(axis="both", labelsize=default_fontsize)
    # ax[1].set_xticks( x_array )
    # ax[1].set_xticklabels( ue_array )
    # ax[1].grid( axis="y", alpha=0.4 )
    # ax[1].legend(ncol=4, loc='upper center', frameon=False, handlelength=1.0, handleheight=1.0, fontsize=default_fontsize - 5)

    plt.tight_layout()
    fig.savefig("nvs-optimal"+ INTRA + FTYPE)

plot_exp_graph()
