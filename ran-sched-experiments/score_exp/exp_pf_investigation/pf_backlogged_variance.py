import re
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt


# Configuration for the number of users per slice
ues_per_slice = [14, 11, 7, 9, 13, 9, 7, 12, 13, 6, 14, 10, 7, 13, 11]

def read_log_file_reverse(filename):
    with open(filename, 'r') as file:
        return file.readlines()[::-1]

def extract_data(lines, start_slice, end_slice):
    pattern = r"app: (\d+) cumu_bytes: (\d+)"
    data = defaultdict(int)
    for line in lines:
        match = re.search(pattern, line)
        if match:
            app_id, cumu_bytes = int(match.group(1)), int(match.group(2))
            slice_index = determine_slice(app_id, ues_per_slice)
            if start_slice <= slice_index <= end_slice:
                if app_id not in data:
                    data[app_id] = cumu_bytes
    return data

def determine_slice(app_id, ues_per_slice):
    total_apps = 0
    for slice_index, num_apps in enumerate(ues_per_slice):
        total_apps += num_apps
        if app_id < total_apps:
            return slice_index
    return -1  # Invalid app_id

def calculate_variance(data, ues_per_slice, start_slice, end_slice):
    variances = []
    for slice_index in range(start_slice, end_slice + 1):
        slice_data = [bytes for app, bytes in data.items() if determine_slice(app, ues_per_slice) == slice_index]
        variances.append(np.var(slice_data))
    return variances

print("Interslice is MT0 *******************************")
filename = "./exp-customize-20slices/max_throughput_0.log"  

lines = read_log_file_reverse(filename)
data = extract_data(lines, 5, 9)
print(min(data.keys()), max(data.keys()))
variances_mt = calculate_variance(data, ues_per_slice, 5, 9)
print(variances_mt)


print("Interslice is PF0 *******************************")
filename = "./exp-customize-20slices/pf_0.log"  
lines = read_log_file_reverse(filename)
data = extract_data(lines, 5, 9)
print(min(data.keys()), max(data.keys()))
variances_pf = calculate_variance(data, ues_per_slice, 5, 9)
print(variances_pf)


print("Interslice is mlwdf0 *******************************")
filename = "./exp-customize-20slices/mlwdf_0.log"  
lines = read_log_file_reverse(filename)
data = extract_data(lines, 5, 9)
print(min(data.keys()), max(data.keys()))
variances_mlwdf = calculate_variance(data, ues_per_slice, 5, 9)
print(variances_mlwdf)


num_slices = 5
# Creating a figure and axis object
fig, ax = plt.subplots()

# Setting the positions and width for the bars
bar_width = 0.25
positions = np.arange(num_slices)

# Plotting the data
bar_mt = ax.bar(positions, variances_mt, bar_width, label='MT')
bar_pf = ax.bar(positions + bar_width, variances_pf, bar_width, label='PF')
bar_mlwdf = ax.bar(positions + 2 * bar_width, variances_mlwdf, bar_width, label='MLWDF')

# Adding labels and title
ax.set_xlabel('Slices')
ax.set_ylabel('Variance of Cumulative Bytes')
ax.set_title('Variance Comparison Across Different Scheduling Algorithms log 0')
ax.set_xticks(positions + bar_width)
ax.set_xticklabels([f'Slice {i+1}' for i in range(num_slices)])
ax.legend()

plt.savefig("pf_backlogged_variance_log0.png")

