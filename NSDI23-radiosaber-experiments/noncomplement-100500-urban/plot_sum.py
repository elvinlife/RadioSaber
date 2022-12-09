#!/usr/bin/python3
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

TIMES=3
INTRA=""
INPUT_DIR="randues"
FTYPE=".pdf"
n_users=450

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
        for line in fin:
            words = line.split(" ")
            if not words[0].isdigit():
                continue
            if int(words[0]) > end_ts:
                break
            if int(words[0]) > begin_ts:
                flow = int(words[2])
                sid = int(words[12])
                cumu_rbs[flow] = int( words[6] ) / (end_ts / 1000 )
                cumu_bytes[flow] = int( words[4] ) / ( end_ts / 1000 )
                flow_to_slice[flow] = sid

    for fid in range(n_users):
        sid = flow_to_slice[fid]
        if sid == -1:
            continue
        slice_cumu_bytes[sid] += cumu_bytes[fid]
        slice_cumu_rbs[sid] += cumu_rbs[fid]

    return slice_cumu_rbs, slice_cumu_bytes

def get_throughput_schemes(dname, n_slices):
    ratio = 8 / (1000 * 1000)
    avg_throughput = {}
    avg_throughput['nvs'] = []
    avg_throughput['maxcell'] = []
    avg_rbs = {}
    avg_rbs['nvs'] = [ 0 for i in range(n_slices) ]
    avg_rbs['maxcell'] = [ 0 for i in range(n_slices) ]

    for i in range(0, TIMES):
        nvs_rbs, nvs_bytes = get_cumubytes( dname + "/nvs_" + INTRA + str(i) + ".log", n_slices )
        maxcell_rbs, maxcell_bytes = get_cumubytes( dname + "/maxcell_" + INTRA + str(i) + ".log", n_slices )
        avg_throughput['nvs'].append( sum(nvs_bytes) * ratio )
        avg_throughput['maxcell'].append( sum(maxcell_bytes) * ratio )
    print(avg_throughput)

    return avg_rbs, avg_throughput

def plot_sum_bandwidth():
    default_font = 24
    n_slices = 20
    _, avg_throughput = get_throughput_schemes( INPUT_DIR, n_slices )
    fig, ax = plt.subplots(figsize=(6,6))
    x_array = np.arange( 0, 2, 1 )
    bw_array = [ np.mean(avg_throughput['nvs']), np.mean(avg_throughput['maxcell']) ]
    bwerr_array = [ np.std(avg_throughput['nvs']), np.std(avg_throughput['maxcell']) ]
    scheme_array = ['NVS', 'RadioSaber']

    ax.grid( axis="y", alpha = 0.4 )
    barlist = ax.bar( x_array, bw_array, width = 0.2 )
    barlist[0].set_color("dimgrey")
    barlist[1].set_color("cornflowerblue")
    ax.errorbar( x_array, bw_array, bwerr_array, fmt=".", capsize=12, color="black" )
    ax.set_xlim( -0.5, 1.5 )
    ax.set_ylim( 0, 300 )

    ax.set_xticks( x_array  )
    ax.set_xticklabels( scheme_array )
    ax.tick_params(axis="y", labelsize=default_font)
    ax.tick_params(axis="x", labelsize=default_font + 2)
    ax.set_ylabel( "Throughput(Mbps)", fontsize=default_font + 4)
    plt.tight_layout()
    fig.savefig("nocomplement_bw" + FTYPE )

plot_sum_bandwidth()
