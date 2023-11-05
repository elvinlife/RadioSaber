#!/usr/bin/python3
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

INTRA=""
TIMES=3
INPUT_DIR="results"
FTYPE=".png"
n_users=204
COLORS=["brown", "dimgrey", "cornflowerblue"]

def get_cdf(data, ratio=0):
    data.sort()
    x_array = []
    y_array = []
    for i, item in enumerate(data):
        if i / len(data) >= ratio:
            x_array.append( item )
            y_array.append( i / len(data) )
    return x_array, y_array

# get the per-second cumulative sent bytes
def get_cumubytes(fname, n_slices):
    begin_ts = 0
    end_ts = 10000
    cumu_bytes = [0 for i in range(n_users)]
    cumu_rbs = [0 for i in range(n_users)]
    flow_to_slice = [-1 for i in range(n_users)]
    slice_cumu_bytes = [0 for i in range(n_slices)]
    slice_cumu_rbs = [0 for i in range(n_slices)]

    with open(fname, "r") as fin:
        # for line in fin:
        #     words = line.split(" ")
        #     if not words[0].isdigit():
        #         continue
        #     if int(words[0]) > end_ts:
        #         break
        #     if int(words[0]) > begin_ts:
        #         flow = int(words[2])
        #         sid = int(words[13])
        #         cumu_rbs[flow] = int( words[7] ) / (end_ts / 1000 )
        #         cumu_bytes[flow] = int( words[5] ) / ( end_ts / 1000 )
        #         flow_to_slice[flow] = sid
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
                flow = int(words[words.index('app:') + 1])
                sid = int(words[words.index('slice:') + 1])
                cumu_rbs[flow] = int(words[words.index('cumu_rbs:') + 1]) / (end_ts / 1000)
                cumu_bytes[flow] = int(words[words.index('cumu_bytes:') + 1]) / (end_ts / 1000)
                flow_to_slice[flow] = sid
                # print(flow, sid, cumu_rbs[flow], cumu_bytes[flow])

    for fid in range(n_users):
        sid = flow_to_slice[fid]
        if sid == -1:
            continue
        slice_cumu_bytes[sid] += cumu_bytes[fid]
        slice_cumu_rbs[sid] += cumu_rbs[fid]
        
    # print(slice_cumu_bytes[sid], " ", sid)
    # print(slice_cumu_rbs[sid], " ",sid)

    return slice_cumu_rbs, slice_cumu_bytes

def get_throughput_perslice(dname, n_slices):
    ratio = 8 / (1000 * 1000)
    avg_throughput = {}
    avg_throughput['mlwdf'] = [0 for i in range(n_slices)]
    avg_throughput['max_throughput'] = [0 for i in range(n_slices)]
    avg_throughput['pf'] = [0 for i in range(n_slices)]
    avg_rbs = {}
    avg_rbs['mlwdf'] = [ 0 for i in range(n_slices) ]
    avg_rbs['max_throughput'] = [ 0 for i in range(n_slices) ]
    avg_rbs['pf'] = [0 for i in range(n_slices) ]

    for i in range(0, TIMES):
        mlwdf_rbs, mlwdf_bytes = get_cumubytes( dname + "/backlogged_m_lwdf_" + str(i) + ".log", n_slices )
        max_throughput_rbs, max_throughput_bytes = get_cumubytes( dname + "/backlogged_max_throughput_"+ str(i) + ".log", n_slices )
        pf_rbs, pf_bytes = get_cumubytes( dname + "/backlogged_proportional_fairness_" + str(i) + ".log", n_slices )
        for j in range(n_slices):
            avg_throughput['mlwdf'][j] += (mlwdf_bytes[j] * ratio)
            avg_throughput['max_throughput'][j] += (max_throughput_bytes[j] * ratio)
            avg_throughput['pf'][j] += (pf_bytes[j] * ratio)
            avg_rbs['mlwdf'][j] += mlwdf_rbs[j]
            avg_rbs['max_throughput'][j] += max_throughput_rbs[j]
            avg_rbs['pf'][j] += pf_rbs[j]

    return avg_rbs, avg_throughput

def get_throughput_total(dname, n_slices):
    ratio = 8 / (1000 * 1000)
    avg_throughput = {}
    avg_throughput['mlwdf'] = []
    avg_throughput['max_throughput'] = []
    avg_throughput['pf'] = []

    for i in range(0, TIMES):
        _, mlwdf_bytes = get_cumubytes( dname + "/backlogged_m_lwdf_" + str(i) + ".log", n_slices )
        _, max_throughput_bytes = get_cumubytes( dname + "/backlogged_max_throughput_"+ str(i) + ".log", n_slices )
        _, pf_bytes = get_cumubytes( dname + "/backlogged_proportional_fairness_" + str(i) + ".log", n_slices )
        avg_throughput['mlwdf'].append( sum(mlwdf_bytes) * ratio )
        avg_throughput['max_throughput'].append( sum(max_throughput_bytes) * ratio )
        avg_throughput['pf'].append( sum(pf_bytes) * ratio )

    return [], avg_throughput

def plot_fairness():
    default_font = 20
    n_slices = 20
    avg_rbs, avg_throughput = get_throughput_perslice( INPUT_DIR, n_slices )
    # print(avg_rbs, avg_throughput)
    # print(len(avg_throughput['max_throughput']))
    # print((avg_throughput['max_throughput']))
    print(avg_rbs["mlwdf"])
    print(avg_rbs["max_throughput"])

    fig, ax = plt.subplots(figsize=(10, 6))
    x_array = np.arange(1, n_slices+1 )
    ax.plot( x_array, avg_throughput['max_throughput'], 'bX--', label="Max_throughput" )
    ax.plot( x_array, avg_throughput['pf'], 'yD--', label="PropF" )
    ax.plot( x_array, avg_throughput['mlwdf'], 'ro--', label="MLWDF" )
    ax.set_xlabel("Slice Id", fontsize=default_font + 4)
    ax.set_ylabel("Throughput(Mbps)", fontsize=default_font + 4)
    ax.set_ylim( bottom = 0, top = 100 )
    ax.set_xticks( [ i for i in range(1, n_slices+1, 2)] )
    ax.tick_params(axis="both", labelsize=default_font)
    ax.grid( axis="y", alpha=0.4 )
    ax.legend(fontsize=default_font + 2, frameon=False, ncol=3, loc="upper center")
    plt.tight_layout()
    fig.savefig("sameweight_tenant_bw" + FTYPE )

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot( x_array, avg_rbs['mlwdf'], 'ro--', label="MLWDF" )
    ax.plot( x_array, avg_rbs['max_throughput'], 'bX--', label="Max_throughput" )
    ax.plot( x_array, avg_rbs['pf'], 'yD--', label="PF" )
    ax.set_xlabel("Slice Id", fontsize=default_font + 4)
    ax.set_ylabel("Resource blocks(per-second)", fontsize=default_font + 4)
    ax.set_ylim( top = 80000, bottom= 60000)
    ax.legend(fontsize=default_font + 2)
    ax.set_xticks( [ i for i in range(1, n_slices+1, 2)] )
    ax.tick_params(axis="both", labelsize=default_font)
    ax.grid( axis="y", alpha=0.4 )
    ax.legend(fontsize=default_font + 2, frameon=False, ncol=3, loc="upper center")
    plt.tight_layout()
    fig.savefig("sameweight_tenant_rbs" + FTYPE )
    #fig.savefig("diffweight_tenant_rbs" + FTYPE )

def plot_sum_bandwidth():
    default_font = 24
    n_slices = 20
    _, avg_throughput = get_throughput_total( INPUT_DIR, n_slices )
    fig, ax = plt.subplots(figsize=(6,6))
    x_array = np.arange( 0, 2, 1 )
    bw_array = [ np.mean(avg_throughput['mlwdf']), np.mean(avg_throughput['max_throughput']) ]
    bwerr_array = [ np.std(avg_throughput['mlwdf']), np.std(avg_throughput['max_throughput']) ]
    scheme_array = ['MLWDF', 'Max_throughput']

    ax.grid( axis="y", alpha = 0.4 )
    barlist = ax.bar( x_array, bw_array, width = 0.2 )
    barlist[0].set_color(COLORS[1])
    barlist[1].set_color(COLORS[2])
    ax.errorbar( x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black" )
    ax.set_xlim( -0.5, 1.5 )
    ax.set_ylim( 0, 700 )

    ax.set_xticks( x_array  )
    ax.set_xticklabels( scheme_array )
    ax.tick_params(axis="y", labelsize=default_font)
    ax.tick_params(axis="x", labelsize=default_font + 2)
    ax.set_ylabel( "Throughput(Mbps)", fontsize=default_font + 4)
    plt.tight_layout()
    fig.savefig("real_bw" + FTYPE )

def plot_together():
    default_font = 20
    n_slices = 20
    avg_rbs, avg_throughput = get_throughput_perslice( INPUT_DIR, n_slices )

    fig, ax = plt.subplots(ncols=3, figsize=(26, 7), gridspec_kw={"width_ratios": [2,2,1]})
    x_array = np.arange(1, n_slices+1 )
    ax[0].plot( x_array, avg_rbs['pf'], 'o--', label="PF", markersize=12, color=COLORS[0] )
    ax[0].plot( x_array, avg_rbs['mlwdf'], 'v--', label="MLWDF", markersize=12, color=COLORS[1] )
    ax[0].plot( x_array, avg_rbs['max_throughput'], '^--', label="Max_throughput", markersize=12, color=COLORS[2] )
    ax[0].set_xlabel("Slice Id", fontsize=default_font + 4)
    ax[0].set_ylabel("Resource blocks(per-second)", fontsize=default_font + 4)
    ax[0].set_ylim( top =  80000, bottom = 75000)
    ax[0].legend(fontsize=default_font + 2)
    ax[0].set_xticks( [ i for i in range(1, n_slices+1, 2)] )
    ax[0].tick_params(axis="both", labelsize=default_font)
    ax[0].grid( axis="y", alpha=0.4 )
    ax[0].legend(fontsize=default_font + 2, frameon=False, ncol=3, loc="upper center")
    ax[0].set_title("(a) per-slice RBs per second", y=-0.3, fontsize=default_font+5)

    ax[1].plot( x_array, avg_throughput['pf'], 'o--', label="PF", markersize=12, color=COLORS[0] )
    ax[1].plot( x_array, avg_throughput['mlwdf'], 'v--', label="MLWDF", markersize=12, color=COLORS[1] )
    ax[1].plot( x_array, avg_throughput['max_throughput'], '^--', label="Max_throughput", markersize=12, color=COLORS[2] )
    ax[1].set_xlabel("Slice Id", fontsize=default_font + 4)
    ax[1].set_ylabel("Throughput(Mbps)", fontsize=default_font + 4)
    ax[1].set_ylim( bottom = 0, top = 100 )
    ax[1].set_xticks( [ i for i in range(1, n_slices+1, 2)] )
    ax[1].tick_params(axis="both", labelsize=default_font)
    ax[1].grid( axis="y", alpha=0.4 )
    ax[1].legend(fontsize=default_font + 2, frameon=False, ncol=3, loc="upper center")
    ax[1].set_title("(b) per-slice throughput", y=-0.3, fontsize=default_font+5)

    _, sum_throughput = get_throughput_total( INPUT_DIR, n_slices )
    x_array = np.arange( 0, 3, 1 )
    bw_array = [ np.mean(sum_throughput['pf']), np.mean(sum_throughput['mlwdf']), np.mean(sum_throughput['max_throughput'])]
    bwerr_array = [ np.std(sum_throughput['pf']), np.std(sum_throughput['mlwdf']), np.std(sum_throughput['max_throughput']) ]
    scheme_array = [ "PF", "MLWDF", "Max_throughput" ]
    barlist = ax[2].bar( x_array, bw_array, width = 0.3 )
    for i in range(3):
        barlist[i].set_color(COLORS[i])
    ax[2].errorbar( x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black" )
    ax[2].set_xlim( -0.5, 2.5 )
    ax[2].set_ylim( 0, 700 )
    ax[2].set_xticks( x_array  )
    ax[2].set_xticklabels( scheme_array )
    ax[2].tick_params(axis="both", labelsize=default_font )
    ax[2].grid( axis="y", alpha = 0.4 )
    # ax[2].set_xlabel( "Schemes", fontsize=default_font + 4)
    ax[2].set_ylabel( "Throughput(Mbps)", fontsize=default_font + 4)
    ax[2].set_title("(c) sum total throughput", y=-0.3, fontsize=default_font+5)

    plt.tight_layout()
    # fig.savefig("diffweight_together" + FTYPE )
    fig.savefig("sameweight_together" + FTYPE )

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plot_fairness()
plot_together()
plot_sum_bandwidth()
