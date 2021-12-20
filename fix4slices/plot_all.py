TIMES=8
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

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
            if int(words[0]) > begin_ts:
                flow = int(words[2])
                cumu_rbs[flow] = int( words[6] )
                cumu_bytes[flow] = int( words[4] )
            if int(words[0]) > end_ts:
                break
    return sum(cumu_bytes)//1000, sum(cumu_rbs)//1000

def get_onecase(dname, ues):
    cumu_bytes = defaultdict(list)
    cumu_rbs = defaultdict(list)
    for i in range(TIMES):
        b, r = get_onetrace( dname + "/nvs_pf" + str(i) + ".log", ues)
        cumu_bytes['nvs'].append(b)
        cumu_rbs['nvs'].append(r)
        b, r = get_onetrace( dname + "/oracle_pf" + str(i) + ".log", ues)
        cumu_bytes['greedy'].append(b)
        cumu_rbs['greedy'].append(r)
        b, r = get_onetrace( dname + "/single_pf" + str(i) + ".log", ues)
        cumu_bytes['single'].append(b)
        cumu_rbs['single'].append(r)
        assert( abs(cumu_rbs['nvs'][-1] - cumu_rbs['greedy'][-1]) < 1 )
        assert( abs(cumu_rbs['nvs'][-1] - cumu_rbs['single'][-1]) < 1 )
    nvs_greedy = []
    greedy_single = []
    for i in range(TIMES):
        if cumu_bytes['nvs'][i] == 0:
            print("%s, %d" % ( dname, i) )
        nvs_greedy.append( cumu_bytes['greedy'][i] / cumu_bytes['nvs'][i] - 1 )
        greedy_single.append( cumu_bytes['greedy'][i] / cumu_bytes['single'][i] - 1)
    return nvs_greedy, greedy_single

def plot_bar():
    dnames = ['2ues', '4ues', '6ues', '8ues']
    x_array = np.array( [2, 4, 6, 8] )
    y1_array = []
    y2_array = []
    yerr1_array = []
    yerr2_array = []
    for i in range(4):
        nvs_greedy, greedy_single = get_onecase( dnames[i], x_array[i]*4)
        y1_array.append( np.mean( nvs_greedy ) )
        yerr1_array.append( np.std( nvs_greedy ) )
        y2_array.append( np.mean( greedy_single ) )
        yerr2_array.append( np.std( greedy_single ) )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar( x_array - 0.25, y1_array, color="cornflowerblue", width=0.5, label="NVS->Greedy")
    ax.bar( x_array + 0.25, y2_array, color="orange", width=0.5, label="Intra->Greedy")
    ax.errorbar(x_array - 0.25, y1_array, yerr1_array, fmt=".", elinewidth=1, capsize=10)
    ax.errorbar(x_array + 0.25, y2_array, yerr2_array, fmt=".", elinewidth=1, capsize=10)
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim( (-0.2, 0.8 ) )
    ax.set_xlabel("# of UEs", fontsize=14)
    ax.set_ylabel("Throughput Gain/Loss", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend()
    fig.savefig("fix10slices.png")

plot_bar()

