import re

def calculate_average_hol_delay(log_filename):
    # Regular expression pattern to match lines in the log file
    line_pattern = re.compile(r'(\d+) app: (\d+) cumu_bytes: (\d+) cumu_rbs: (\d+) hol_delay: ([\d.]+) user: (\d+) slice: (\d+)')

    # Initialize a dictionary to store the total hol_delay and packet count for each slice
    slice_data = {}

    with open(log_filename, 'r') as file:
        for line in file:
            match = line_pattern.match(line.strip())
            if match:
                hol_delay = float(match.group(5))
                slice_id = int(match.group(7))

                if slice_id not in slice_data:
                    slice_data[slice_id] = {'total_hol_delay': 0, 'packet_count': 0}

                slice_data[slice_id]['total_hol_delay'] += hol_delay
                slice_data[slice_id]['packet_count'] += 1

    # Calculate and print average hol_delay for each slice
    average_hol_delays = {}
    for slice_id, data in slice_data.items():
        if slice_id in range(10, 15):
            average_hol_delays[slice_id] = data['total_hol_delay'] / data['packet_count']
            print(f"Slice {slice_id}: Average HOL Delay = {average_hol_delays[slice_id]}")

    return average_hol_delays

# Example usage
print("PF")
log_filename = './exp-customize-20slices/pf_1.log'
calculate_average_hol_delay(log_filename)

print("MT")
log_filename = './exp-customize-20slices/max_throughput_1.log'
calculate_average_hol_delay(log_filename)

print("mlwdf")
log_filename = './exp-customize-20slices/mlwdf_1.log'
calculate_average_hol_delay(log_filename)


