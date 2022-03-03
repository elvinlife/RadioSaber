TIMES=3
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

ues_per_slice = 4

def get_JIndex_array(fname, n_users):
    begin_ts = 9000
    end_ts = 10000
    # cumu_bytes = [0 for i in range(n_users)]
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
                # cumu_bytes[flow] = int( words[4] )
    # slice_rbs = [0 for i in range(n_users//ues_per_slice)]
    rbs_2d = [cumu_rbs[i:i+ues_per_slice] for i in range(0, len(cumu_rbs), ues_per_slice)]
    # print(rbs_2d)
    J_index_array = []
    for i in range(len(rbs_2d)):
        J_index_array.append(calculate_JIndex(rbs_2d[i]))
    # print(J_index_array)
    return J_index_array

def get_inter_JIndex(fname, n_users):
    begin_ts = 9000
    end_ts = 10000
    # cumu_bytes = [0 for i in range(n_users)]
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
                # cumu_bytes[flow] = int( words[4] )
    slice_rbs = [0 for _ in range(n_users//ues_per_slice)]
    for i in range(n_users//ues_per_slice):
        for j in range(ues_per_slice):
            slice_rbs[i] += cumu_rbs[i*ues_per_slice+j]
    return calculate_JIndex(slice_rbs)


def calculate_JIndex(slice_rbs):
    numberOfSlices = len(slice_rbs)
    sum_ = sum(slice_rbs)
    den = sum([i**2 for i in slice_rbs])
    den = den * numberOfSlices
    return sum_**2/den



def get_JainsIndex(dname, ues):
    nvs_index, greedy_index, subopt_index = [], [], []
    nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = [], [], [], [], [], [] # number of slices
    for i in range(TIMES):
        b_nvs_J_index = get_JIndex_array(dname + "/nvs_pf" + str(i) + ".log", ues)
        nvs_index.append(b_nvs_J_index)

        b_greedy_J_index = get_JIndex_array(dname + "/greedy_pf" + str(i) + ".log", ues)
        greedy_index.append(b_greedy_J_index)

        b_suboptimal_J_index = get_JIndex_array(dname + "/subopt_pf" + str(i) + ".log", ues)
        subopt_index.append(b_suboptimal_J_index)
    # J_index["nvs"] = J_index["nvs"]/TIMES
    # J_index["greedy"] = J_index["greedy"]/TIMES
    # J_index["subopt"] = J_index["subopt"]/TIMES

    for i in range(len(nvs_index[0])):
        nvs_temp = []
        greedy_temp = []
        subopt_temp = []
        for j in range(TIMES):
            nvs_temp.append(nvs_index[j][i])
            greedy_temp.append(greedy_index[j][i])
            subopt_temp.append(subopt_index[j][i])
        nvs_mean.append(np.mean(nvs_temp))
        nvs_err.append(np.std(nvs_temp))
        greedy_mean.append(np.mean(greedy_temp))
        greedy_err.append(np.std(greedy_temp))
        subopt_mean.append(np.mean(subopt_temp))
        subopt_err.append(np.std(subopt_temp))

    return nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err


def get_inter_JainsIndex(dname, ues):
    nvs_index, greedy_index, subopt_index = 0,0,0
    # nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = [], [], [], [], [], [] # number of slices
    for i in range(TIMES):
        b_nvs_J_index = get_inter_JIndex(dname + "/nvs_pf" + str(i) + ".log", ues)
        nvs_index+=b_nvs_J_index

        b_greedy_J_index = get_inter_JIndex(dname + "/greedy_pf" + str(i) + ".log", ues)
        greedy_index+=b_greedy_J_index

        b_suboptimal_J_index = get_inter_JIndex(dname + "/subopt_pf" + str(i) + ".log", ues)
        subopt_index+=b_suboptimal_J_index
    nvs_index/=TIMES
    greedy_index/=TIMES
    subopt_index/=TIMES
    return nvs_index, greedy_index, subopt_index




def plot_bar_throughput():
    #dnames = ['4slices', '10slices', '20slices', '40slices', '100slices']
    #x_array = np.arange(5)
    #slice_array = [4, 10, 20, 40, 100]

    dnames = ['4slices', '10slices', '20slices', '40slices']
    x_array = np.arange(8, step=2)
    slice_array = [4, 10, 20, 40]
    y1_array = []
    y2_array = []
    y3_array = []
    y4_array = []
    y5_array = []
    yerr1_array = []
    yerr2_array = []
    yerr3_array = []
    yerr4_array = []
    yerr5_array = []
    for i in range(len(dnames)):
        cumu_bytes = get_throughput( dnames[i], slice_array[i]*4 )
        y1_array.append( np.mean( cumu_bytes['nvs'] ) )
        yerr1_array.append( np.std( cumu_bytes['nvs'] ) )
        y2_array.append( np.mean( cumu_bytes['greedy'] ) )
        yerr2_array.append( np.std( cumu_bytes['greedy'] ) )
        y3_array.append( np.mean( cumu_bytes['subopt'] ) )
        yerr3_array.append( np.std( cumu_bytes['subopt'] ) )
        y4_array.append( np.mean( cumu_bytes['maxcell'] ) )
        yerr4_array.append( np.std( cumu_bytes['maxcell'] ) )
        y5_array.append( np.mean( cumu_bytes['vogel'] ) )
        yerr5_array.append( np.std( cumu_bytes['vogel'] ) )
    fig, ax = plt.subplots(figsize=(8, 6))

    #ax.bar( x_array - 0.3, y1_array, color="cornflowerblue", width=0.3, label="NVS")
    #ax.bar( x_array, y2_array, color="orange", width=0.3, label="Greedy")
    #ax.bar( x_array + 0.3, y3_array, color="green", width=0.3, label="Subopt")
    #ax.errorbar(x_array - 0.3, y1_array, yerr1_array, fmt=".", elinewidth=0.3, capsize=10)
    #ax.errorbar(x_array, y2_array, yerr2_array, fmt=".", elinewidth=0.3, capsize=10)
    #ax.errorbar(x_array + 0.3, y3_array, yerr3_array, fmt=".", elinewidth=0.3, capsize=10)

    ax.bar( x_array - 0.6, y1_array, color="cornflowerblue", width=0.3, label="NVS")
    ax.bar( x_array - 0.3, y2_array, color="orange", width=0.3, label="Greedy")
    ax.bar( x_array + 0.0, y3_array, color="green", width=0.3, label="Subopt")
    ax.bar( x_array + 0.3, y4_array, color="black", width=0.3, label="MaxCell")
    ax.bar( x_array + 0.6, y5_array, color="blue", width=0.3, label="Vogel")
    ax.errorbar(x_array - 0.6, y1_array, yerr1_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array - 0.3, y2_array, yerr2_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array - 0.0, y3_array, yerr3_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array + 0.3, y4_array, yerr4_array, fmt=".", elinewidth=0.3, capsize=10)
    ax.errorbar(x_array + 0.6, y5_array, yerr5_array, fmt=".", elinewidth=0.3, capsize=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim( bottom = 1 )
    ax.set_xlabel("# of slices", fontsize=14)
    ax.set_ylabel("Throughput(Mbps)", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xticks( x_array )
    ax.set_xticklabels( slice_array )
    ax.legend(loc="lower left")
    fig.savefig("fix4ues.png")

# plot_bar_throughput()

def plot_JainIndex():
    fig, axs = plt.subplots(2,2)
    fig.suptitle("fixing 4 ues")
    # x1 = 
    # fig.savefig("JainIndexFixing4ues.png")


    dnames = ['4slices', '10slices', '20slices', '40slices']
    slice_array = [4, 10, 20, 40]
    
    x1_array = np.arange(4)
    x2_array = np.arange(10)
    x3_array = np.arange(20)
    x4_array = np.arange(40)

    nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = get_JainsIndex(dnames[0], slice_array[0]*4)
    axs[0,0].bar(x1_array, nvs_mean, color="cornflowerblue", width=0.1, label="NVS")
    axs[0,0].errorbar(x1_array, nvs_mean, nvs_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,0].bar(x1_array+0.1, greedy_mean, color="orange", width=0.1, label="Greedy")
    axs[0,0].errorbar(x1_array+0.1, greedy_mean, greedy_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,0].bar(x1_array+0.2, subopt_mean, color="green", width=0.1, label="Subopt")
    axs[0,0].errorbar(x1_array+0.2, subopt_mean, subopt_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,0].title.set_text("4slices")
    axs[0,0].legend(loc="lower left")

    nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = get_JainsIndex(dnames[1], slice_array[1]*4)
    axs[0,1].bar(x2_array, nvs_mean, color="cornflowerblue", width=0.1, label="NVS")
    axs[0,1].errorbar(x2_array, nvs_mean, nvs_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,1].bar(x2_array+0.1, greedy_mean, color="orange", width=0.1, label="Greedy")
    axs[0,1].errorbar(x2_array+0.1, greedy_mean, greedy_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,1].bar(x2_array+0.2, subopt_mean, color="green", width=0.1, label="Subopt")
    axs[0,1].errorbar(x2_array+0.2, subopt_mean, subopt_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[0,1].title.set_text("10slices")

    nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = get_JainsIndex(dnames[2], slice_array[2]*4)
    axs[1,0].bar(x3_array, nvs_mean, color="cornflowerblue", width=0.1, label="NVS")
    axs[1,0].errorbar(x3_array, nvs_mean, nvs_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,0].bar(x3_array+0.1, greedy_mean, color="orange", width=0.1, label="Greedy")
    axs[1,0].errorbar(x3_array+0.1, greedy_mean, greedy_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,0].bar(x3_array+0.2, subopt_mean, color="green", width=0.1, label="Subopt")
    axs[1,0].errorbar(x3_array+0.2, subopt_mean, subopt_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,0].title.set_text("20slices")

    nvs_mean, nvs_err, greedy_mean, greedy_err, subopt_mean, subopt_err = get_JainsIndex(dnames[3], slice_array[3]*4)
    axs[1,1].bar(x4_array, nvs_mean, color="cornflowerblue", width=0.1, label="NVS")
    axs[1,1].errorbar(x4_array, nvs_mean, nvs_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,1].bar(x4_array+0.1, greedy_mean, color="orange", width=0.1, label="Greedy")
    axs[1,1].errorbar(x4_array+0.1, greedy_mean, greedy_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,1].bar(x4_array+0.2, subopt_mean, color="green", width=0.1, label="Subopt")
    axs[1,1].errorbar(x4_array+0.2, subopt_mean, subopt_err, fmt=".", elinewidth=0.1, capsize=1)
    axs[1,1].title.set_text("40slices")

    fig.savefig("JainIndexFixing4ues.png")


    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle("fixing 4 ues InterSlice Jain Index")
    x_array = [i for i in range(4)]
    y_array_nvs = []
    y_array_greedy = []
    y_array_subopt = []
    for i in range(4):
        nvs_index, greedy_index, subopt_index = get_inter_JainsIndex(dnames[i], slice_array[i]*4)
        y_array_nvs.append(nvs_index)
        y_array_greedy.append(greedy_index)
        y_array_subopt.append(subopt_index)
    ax.plot(x_array, y_array_nvs, label = "NVS")
    ax.plot(x_array, y_array_greedy, label = "Greedy")
    ax.plot(x_array, y_array_subopt, label = "Subopt")
    ax.legend()

    fig.savefig("JainIndexFixing4uesInter.png")




plot_JainIndex()
