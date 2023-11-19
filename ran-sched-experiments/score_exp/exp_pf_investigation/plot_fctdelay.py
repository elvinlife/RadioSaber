#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from collections import defaultdict
import numpy as np

BEGIN_IDX=0
INTRA=""
INPUT_DIR="exp-customize-20slices"
FTYPE=".png"
COLORS=["brown", "dimgrey", "cornflowerblue"]

def get_cdf(data, ratio=0, upperbound=100000):
    data.sort()
    x_array = []
    y_array = []
    for i, item in enumerate(data):
        if i / len(data) >= ratio:
            if item > upperbound:
                break
            x_array.append( item )
            y_array.append( i / len(data) )
    return x_array, y_array

def get_hol(fname, slice_begin, slice_end):
    hol_array = []
    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if not words[0].isdigit():
                continue
            flow = int(words[2])
            sid = int(words[12])
            if sid >= slice_begin and sid <= slice_end:
                hol_array.append( float(words[8] ) )
    return hol_array

def get_fct(fname, slice_begin, slice_end, priority_only):
    flow_to_slice = {}
    fct_captured = {}
    # only analyze flows generated before ts_shoot
    ts_shoot = 10000
    is_after_ts_shoot = False
    fct_array = []
    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if words[0].isdigit():
                app_id = int(words[2])
                flow_to_slice[app_id] = int(words[12])
    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if words[0].isdigit():
                if int(words[0]) > ts_shoot:
                    is_after_ts_shoot = True
            if words[0] == "ipflow":
                # even number is higher priority
                app_id = int(words[3])
                # skip flows not belong to this slice
                if app_id in flow_to_slice:
                    if flow_to_slice[app_id] < slice_begin or flow_to_slice[app_id] > slice_end:
                        continue
                else:
                    print("%s, flow: %d never scheduled" % (fname, app_id) )

                if words[1] == "start" and not is_after_ts_shoot:
                    key =( app_id, int(words[5]) )
                    fct_captured[key] = -1
                if words[1] == "end":
                    # skip non-priority flows
                    if priority_only and int(words[11]) == 0:
                        continue
                    key =( app_id, int(words[5]) )
                    if key in fct_captured:
                        fct_captured[key] = float( words[7] )
    for k, v in fct_captured.items():
        if v != -1:
            fct_array.append( v )
    print( "%s: finished %d flows before %dms" % ( fname, len(fct_array), ts_shoot) )
    return fct_array

def get_throughput(fname, slice_begin, slice_end):
    perflow_throughput = {}
    cumu_throughput = {}
    begin_ts = 20000
    end_ts = 22000
    for i in range(slice_begin, slice_end):
        perflow_throughput[i] = {}
        cumu_throughput[i] = 0
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
                if sid >= slice_begin and sid < slice_end:
                    perflow_throughput[sid][flow] = int( words[4] ) / (end_ts / 1000 ) * 8 / (1000 * 1000)
    for i in range(slice_begin, slice_end):
        for k in perflow_throughput[i].keys():
            cumu_throughput[i] += perflow_throughput[i][k]
    return [v for v in cumu_throughput.values()]

def get_fct_schemes(dname, slice_begin, slice_end, priority_only):
    all_fct = {}
    for i in range(BEGIN_IDX, BEGIN_IDX+1):
        all_fct['mt'] = get_fct( dname + "/max_throughput_" + str(i) + ".log", slice_begin, slice_end, priority_only )
        all_fct['mlwdf'] = get_fct( dname + "/mlwdf_" + str(i) + ".log", slice_begin, slice_end, priority_only )
        all_fct['pf'] = get_fct( dname + "/pf_" + INTRA + str(i) + ".log", slice_begin, slice_end, priority_only )
    return all_fct

def get_hol_schemes(dname, slice_begin, slice_end):
    all_hol = {}
    for i in range(BEGIN_IDX, BEGIN_IDX+1):
        all_hol['mt'] = get_hol( dname + "/max_throughput_" + str(i) + ".log", slice_begin, slice_end)
        all_hol['mlwdf'] = get_hol( dname + "/mlwdf_" + str(i) + ".log", slice_begin, slice_end )
        all_hol['pf'] = get_hol( dname + "/pf_" + INTRA + str(i) + ".log", slice_begin, slice_end )
    return all_hol

def get_throughput_schemes(dname, slice_begin, slice_end):
    all_throughput = {}
    for i in range(BEGIN_IDX, BEGIN_IDX+1):
        all_throughput['mt'] = np.mean(
                get_throughput( dname + "/max_throughput_" + str(i) + ".log", slice_begin, slice_end ) )
        all_throughput['mlwdf'] = np.mean(
                get_throughput( dname + "/mlwdf_" + str(i) + ".log", slice_begin, slice_end ) )
        all_throughput['pf'] = np.mean(
                get_throughput( dname + "/pf_" + INTRA + str(i) + ".log", slice_begin, slice_end ) )
    return all_throughput

def print_throughput(slice_begin, slice_end):
    all_throughput = get_throughput_schemes(
            INPUT_DIR, slice_begin, slice_end )
    print(all_throughput)


def plot_hol_delay(ofname, slice_begin, slice_end):
    default_font = 20
    fig, ax = plt.subplots(figsize=(9, 5))
    all_hol = get_hol_schemes( INPUT_DIR, slice_begin, slice_end )
    mt_x, mt_y = get_cdf( all_hol['mt'], 0.0, 5 )
    mlwdf_x, mlwdf_y = get_cdf( all_hol['mlwdf'], 0.0, 5 )
    pf_x, pf_y = get_cdf( all_hol['pf'], 0.0, 5 )
    
    print(np.mean(all_hol['mt']), np.mean(all_hol['mlwdf']), np.mean(all_hol['pf']))

    ax.plot( pf_x, pf_y, "y--", label="PF", color=COLORS[0], linewidth=2.5 )
    ax.plot( mt_x, mt_y, "r--", label="MT", color=COLORS[1], linewidth=2.5 )
    ax.plot( mlwdf_x, mlwdf_y, "b--", label="MLWDF", color=COLORS[2], linewidth=2.5 )
    ax.set_xlabel("Queueing Delay(s)", fontsize=default_font + 4)
    ax.set_ylabel("Ratio", fontsize=default_font + 4)
    ax.legend(fontsize=default_font + 2)
    ax.tick_params(axis="both", labelsize=default_font)
    ax.grid( axis="both", alpha=0.2 )
    plt.tight_layout()
    fig.savefig( ofname + FTYPE )

print_throughput( 0, 4 )
plot_hol_delay( "cdf-hol-delay", 10, 14 )
