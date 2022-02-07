TIMES=1
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def get_cdf(data, ratio=1):
    data.sort()
    x_array = []
    y_array = []
    for i, item in enumerate(data):
        x_array.append( i / len(data) )
        y_array.append( item )
    return x_array, y_array

def get_onetrace(fname, n_users):
    cumu_delays = [0 for i in range(n_users)]
    delay_counts = [0 for i in range(n_users)]
    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if not words[0].isdigit():
                continue
            flow = int(words[2])
            cumu_delays[flow] += float(words[8])
            delay_counts[flow] += 1
    avg_delays = []
    for i in range(n_users):
        avg_delays.append( cumu_delays[i] / delay_counts[i] * 1000 )
    return avg_delays

def get_avg_delay(dname, ues):
    avg_delays = defaultdict(list)
    for i in range(TIMES):
        delays = get_onetrace( dname + "/nvs_pf" + str(i) + ".log", ues)
        avg_delays['nvs'].append( delays )
        delays = get_onetrace( dname + "/greedy_pf" + str(i) + ".log", ues)
        avg_delays['greedy'].append( delays )
    return avg_delays

def plot_hol_cdf():
    dnames = ["50slices", "25slices", "10slices"]
    all_ues = 100
    for i in range(len(dnames)):
        fig, ax = plt.subplots(figsize=(8, 6))
        avg_delays = get_avg_delay(dnames[i], all_ues)
        nvs_x, nvs_y = get_cdf( avg_delays['nvs'][0] )
        greed_x, greed_y = get_cdf( avg_delays['greedy'][0] )
        ax.plot( nvs_x, nvs_y, "b--", label="NVS")
        ax.plot( greed_x, greed_y, "r--", label="Greedy")
        ax.set_xlabel("Ratio")
        ax.set_ylabel("HOL Delay(ms)")
        ax.legend()
        fig.savefig( dnames[i] + ".png" )

plot_hol_cdf()
