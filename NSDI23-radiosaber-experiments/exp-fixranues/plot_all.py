#!/usr/bin/python3
TIMES=1
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
import numpy as np
INTRA="mt"
LABEL="MT"
FTYPE=".pdf"
COLORS=["dimgrey", "brown", "cornflowerblue", "wheat"]

def get_onetrace(fname, n_users):
    begin_ts = 0
    end_ts = 8000
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
    for i in range(TIMES):
        b, _ = get_onetrace( dname + "/nvs_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['nvs'].append(b * ratio)
        b, _ = get_onetrace( dname + "/maxcell_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['maxcell'].append(b * ratio)
        b, _ = get_onetrace( dname + "/sequential_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['sequential'].append(b * ratio)
        b, _ = get_onetrace( dname + "/upperbound_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['upperbound'].append(b * ratio)
    return cumu_bytes

def plot_exp1_graph():
    default_fontsize = 24
    # dnames = ['5slices', '10slices', '15slices', '20slices']
    dnames = ['5slices-ip', '10slices-ip', '15slices-ip', '20slices-ip']
    x_array = np.arange(0, 4, 1)
    slice_array = [5, 10, 15, 20]
    y1_array = []
    y3_array = []
    yerr1_array = []
    yerr3_array = []
    for i in range(len(dnames)):
        cumu_bytes = get_throughput( dnames[i], slice_array[i]*15)
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y3_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr3_array.append( np.std( cumu_bytes['maxcell'] ) )
    fig, ax = plt.subplots(figsize=(8, 5))
    bar_width = 0.2
    ax.bar( x_array - 0.5 * bar_width, y1_array, width=bar_width, label="NVS", color=COLORS[0])
    ax.bar( x_array + 0.5 * bar_width, y3_array, width=bar_width, label="RadioSaber", color=COLORS[2])
    ax.errorbar(x_array - 0.5 * bar_width, y1_array, yerr1_array, fmt=".", capsize=10, color="black")
    ax.errorbar(x_array + 0.5 * bar_width, y3_array, yerr3_array, fmt=".", capsize=10, color="black")

    for i in range(4):
        print("%f" % (y3_array[i] / y1_array[i] - 1) )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 300)
    ax.set_xlabel("Num of Slices", fontsize=default_fontsize + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    ax.tick_params(axis="both", labelsize=default_fontsize - 2)
    ax.grid( alpha=0.4, axis="y" )
    ax.set_xticks( x_array )
    ax.set_xticklabels( slice_array )
    ax.legend(ncol=2, loc='upper center', frameon=False,\
            handlelength=1.0, handleheight=1.0,\
            bbox_to_anchor=(0.5, 1.05),\
            fontsize=default_fontsize + 2)
    plt.tight_layout()

    # fig.savefig("exp1-fixrandues-" + INTRA + FTYPE)
    fig.savefig("exp1-fixrandues-ip-" + INTRA + FTYPE)

def plot_exp2_graph():
    default_fontsize = 22
    dnames = ['5slices', '10slices', '15slices', '20slices']
    #dnames = ['5slices-ip', '10slices-ip', '15slices-ip', '20slices-ip']
    x_array = np.arange(0, 4, 1)
    slice_array = [5, 10, 15, 20]
    y1_array = []
    y2_array = []
    y3_array = []
    y4_array = []
    yerr1_array = []
    yerr2_array = []
    yerr3_array = []
    yerr4_array = []
    for i in range(len(dnames)):
        cumu_bytes = get_throughput( dnames[i], slice_array[i]*15)
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['sequential'] ) )
        yerr2_array.append( np.std( cumu_bytes['sequential'] ) )
        y3_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr3_array.append( np.std( cumu_bytes['maxcell'] ) )
        y4_array.append( np.mean( cumu_bytes['upperbound'] ) )
        yerr4_array.append( np.std( cumu_bytes['upperbound'] ) )
    fig, ax = plt.subplots(figsize=(13, 5))
    bar_width = 0.16
    ax.bar( x_array - 1.5 * bar_width, y1_array, width=bar_width, label="NVS", color=COLORS[0], hatch="", alpha=0.99 )
    ax.bar( x_array - 0.5 * bar_width, y2_array, width=bar_width, label="Sequential", color=COLORS[1], hatch="", alpha=0.99 )
    ax.bar( x_array + 0.5 * bar_width, y3_array, width=bar_width, label="RadioSaber", color=COLORS[2], hatch="", alpha=0.99 )
    ax.bar( x_array + 1.5 * bar_width, y4_array, width=bar_width, label="Upperbound", color=COLORS[3], hatch="", alpha=0.99 )
    ax.errorbar(x_array - 1.5 * bar_width, y1_array, yerr1_array, fmt=".",  capsize=12, color="black")
    ax.errorbar(x_array - 0.5 * bar_width, y2_array, yerr2_array, fmt=".",  capsize=12, color="black")
    ax.errorbar(x_array + 0.5 * bar_width, y3_array, yerr3_array, fmt=".",  capsize=12, color="black")
    ax.errorbar(x_array + 1.5 * bar_width, y4_array, yerr4_array, fmt=".",  capsize=12, color="black")

    for i in range(4):
        print("%d slices: radiosaber %f higher than sequential; %f lower than upperbound" %  \
                ( slice_array[i], y3_array[i] / y2_array[i] - 1, \
                    1 - y3_array[i] / y4_array[i] ) )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(100, 400)
    ax.set_xlabel("Num of Slices", fontsize=default_fontsize + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_fontsize + 4)
    ax.tick_params(axis="both", labelsize=default_fontsize)
    ax.grid( axis="y", alpha=0.4 )
    ax.set_xticks( x_array )
    ax.set_xticklabels( slice_array )
    ax.legend(ncol=4, loc='upper center', frameon=False, handlelength=1.0, handleheight=1.0, fontsize=default_fontsize )
    plt.tight_layout()

    #fig.savefig("exp2-fixrandues-ip-" + INTRA + FTYPE)
    fig.savefig("exp2-fixrandues-" + INTRA + FTYPE)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#plot_exp1_graph()
plot_exp2_graph()
