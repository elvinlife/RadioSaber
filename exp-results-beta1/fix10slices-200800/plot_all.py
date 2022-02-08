#!/usr/bin/python3
TIMES=3
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

INTRA = "mt"

def get_onetrace(fname, n_users):
    begin_ts = 9000
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
                cumu_rbs[flow] = int( words[6] )
                cumu_bytes[flow] = int( words[4] )
    return sum(cumu_bytes)//1000, sum(cumu_rbs)//1000

def get_ratio(dname, ues):
    cumu_bytes = defaultdict(list)
    cumu_rbs = defaultdict(list)
    for i in range(TIMES):
        b, r = get_onetrace( dname + "/nvs_pf" + str(i) + ".log", ues)
        cumu_bytes['nvs'].append(b)
        cumu_rbs['nvs'].append(r)
        b, r = get_onetrace( dname + "/greedy_pf" + str(i) + ".log", ues)
        cumu_bytes['greedy'].append(b)
        cumu_rbs['greedy'].append(r)
        b, r = get_onetrace( dname + "/maxcell_pf" + str(i) + ".log", ues)
        cumu_bytes['maxcell'].append(b)
        cumu_rbs['maxcell'].append(r)
        if abs(cumu_rbs['nvs'][-1] - cumu_rbs['maxcell'][-1]) > 1:
            print("%s, %d, %d, %d" % ( dname, cumu_rbs['nvs'][-1], cumu_rbs['maxcell'][-1], i))
        assert( abs(cumu_rbs['nvs'][-1] - cumu_rbs['greedy'][-1]) <= 1 )
        assert( abs(cumu_rbs['nvs'][-1] - cumu_rbs['maxcell'][-1]) <= 1 )
    nvs_greedy = []
    greedy_maxcell = []
    for i in range(TIMES):
        nvs_greedy.append( cumu_bytes['greedy'][i] / cumu_bytes['nvs'][i] )
        greedy_maxcell.append( cumu_bytes['maxcell'][i] / cumu_bytes['nvs'][i] )
    return nvs_greedy, greedy_maxcell

def get_throughput(dname, ues):
    ratio = 12* 1000 / 8
    cumu_bytes = defaultdict(list)
    for i in range(TIMES):
        b, _ = get_onetrace( dname + "/nvs_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['nvs'].append(b/ratio)
        b, _ = get_onetrace( dname + "/greedy_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['greedy'].append(b/ratio)
        b, _ = get_onetrace( dname + "/subopt_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['subopt'].append(b/ratio)
        b, _ = get_onetrace( dname + "/upperbound_" + INTRA + str(i) + ".log", ues)
        cumu_bytes['single'].append(b/ratio)
    return cumu_bytes

def plot_bar_ratio():
    dnames = ['2slices', '4slices', '6slices', '8slices', '10slices']
    x_array = np.array( [2, 4, 6, 8, 10] )
    y1_array = []
    y2_array = []
    yerr1_array = []
    yerr2_array = []
    for i in range(len(dnames)):
        nvs_greedy, greedy_maxcell = get_ratio( dnames[i], x_array[i]*4 )
        y1_array.append( np.mean( nvs_greedy ) )
        yerr1_array.append( np.std( nvs_greedy ) )
        y2_array.append( np.mean( greedy_maxcell ) )
        yerr2_array.append( np.std( greedy_maxcell ) )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar( x_array - 0.25, y1_array, color="cornflowerblue", width=0.5, label="NVS->Greedy")
    ax.bar( x_array + 0.25, y2_array, color="orange", width=0.5, label="NVS->Optimal")
    ax.errorbar(x_array - 0.25, y1_array, yerr1_array, fmt=".", elinewidth=1, capsize=10)
    ax.errorbar(x_array + 0.25, y2_array, yerr2_array, fmt=".", elinewidth=1, capsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim( bottom = 1 )
    ax.set_xlabel("# of Slices", fontsize=14)
    ax.set_ylabel("Throughput Gain/Loss", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend()
    fig.savefig("fix4ues.png")

#def plot_bar_throughput():
#    dnames = ['2ues', '4ues', '6ues', '8ues', '10ues']
#    x_array = np.array( [2, 4, 6, 8, 10] )
#    y1_array = []
#    y2_array = []
#    y3_array = []
#    y4_array = []
#    y5_array = []
#    yerr1_array = []
#    yerr2_array = []
#    yerr3_array = []
#    yerr4_array = []
#    yerr5_array = []
#    for i in range(len(dnames)):
#        cumu_bytes = get_throughput( dnames[i], x_array[i]*10 )
#        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
#        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
#        y2_array.append( np.mean( cumu_bytes['greedy'] ) )
#        yerr2_array.append( np.std( cumu_bytes['greedy'] ) )
#        y3_array.append( np.mean( cumu_bytes['subopt'] ) )
#        yerr3_array.append( np.std( cumu_bytes['subopt'] ) )
#        y4_array.append( np.mean( cumu_bytes['maxcell'] ) )
#        yerr4_array.append( np.std( cumu_bytes['maxcell'] ) )
#        y5_array.append( np.mean( cumu_bytes['vogel'] ) )
#        yerr5_array.append( np.std( cumu_bytes['vogel'] ) )
#    fig, ax = plt.subplots(figsize=(8, 6))
#    ax.bar( x_array - 0.6, y1_array, color="cornflowerblue", width=0.3, label="NVS")
#    ax.bar( x_array - 0.3, y2_array, color="orange", width=0.3, label="Greedy")
#    ax.bar( x_array , y3_array, color="green", width=0.3, label="Subopt")
#    ax.bar( x_array + 0.3, y4_array, color="yellow", width=0.3, label="MaxCell")
#    ax.bar( x_array + 0.6, y5_array, color="red", width=0.3, label="Vogel")
#    ax.errorbar(x_array - 0.6, y1_array, yerr1_array, fmt=".", elinewidth=0.3, capsize=10)
#    ax.errorbar(x_array - 0.3, y2_array, yerr2_array, fmt=".", elinewidth=0.3, capsize=10)
#    ax.errorbar(x_array , y3_array, yerr3_array, fmt=".", elinewidth=0.3, capsize=10)
#    ax.errorbar(x_array + 0.3, y4_array, yerr4_array, fmt=".", elinewidth=0.3, capsize=10)
#    ax.errorbar(x_array + 0.6, y5_array, yerr5_array, fmt=".", elinewidth=0.3, capsize=10)
#    ax.spines["top"].set_visible(False)
#    ax.spines["right"].set_visible(False)
#    ax.set_ylim( bottom = 1 )
#    ax.set_xlabel("# of ues", fontsize=14)
#    ax.set_ylabel("Throughput(Mbps)", fontsize=14)
#    ax.tick_params(axis="both", labelsize=12)
#    ax.legend(loc="lower left")
#    fig.savefig("fix10slices.png")

def plot_bar_throughput():
    dnames = ['4ues', '10ues', '20ues', '40ues']
    x_array = np.arange(0, 8, 2)
    ue_array = [4, 10, 20, 40]
    y1_array = []
    y2_array = []
    y3_array = []
    y4_array = []
    yerr1_array = []
    yerr2_array = []
    yerr3_array = []
    yerr4_array = []
    for i in range(len(dnames)):
        cumu_bytes = get_throughput( dnames[i], ue_array[i]*10 )
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['greedy'] ) )
        yerr2_array.append( np.std( cumu_bytes['greedy'] ) )
        y3_array.append( np.mean( cumu_bytes['subopt'] ) )
        yerr3_array.append( np.std( cumu_bytes['subopt'] ) )
        y4_array.append( np.mean( cumu_bytes['single'] ) )
        yerr4_array.append( np.std( cumu_bytes['single'] ) )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar( x_array - 0.6, y1_array, color="cornflowerblue", width=0.4, label="NVS")
    ax.bar( x_array - 0.2, y2_array, color="orange", width=0.4, label="Greedy")
    ax.bar( x_array + 0.2, y3_array, color="green", width=0.4, label="Subopt")
    ax.bar( x_array + 0.6, y4_array, color="blue", width=0.4, label="UpperBound")
    ax.errorbar(x_array - 0.6, y1_array, yerr1_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array - 0.2, y2_array, yerr2_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array + 0.2, y3_array, yerr3_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array + 0.6, y4_array, yerr4_array, fmt=".", elinewidth=0.3, capsize=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim( bottom = 1 )
    ax.set_xlabel("# of ues", fontsize=14)
    ax.set_ylabel("Throughput(Mbps)", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xticks( x_array )
    ax.set_xticklabels( ue_array )
    ax.legend(loc="lower left")
    fig.savefig(INTRA + ".png")

plot_bar_throughput()
