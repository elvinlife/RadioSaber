import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Replace with the actual path where your log files are located.
log_files_path = './exp-customize-20slices/*.log'

# This function reads the last 20 lines of a file and extracts the scores.
def extract_last_20_scores(file_path):
    scores = []
    with open(file_path, 'rb') as file:
        # Go to the end of the file
        file.seek(0, os.SEEK_END)
        
        # We are at the end of the file, start moving backwards
        end_byte = file.tell()
        lines_to_go = 15
        buffer = bytearray()
        while lines_to_go > 0 and file.tell() > 0:
            file.seek(-1, os.SEEK_CUR)
            char = file.read(1)
            # When we find a line break, we process the current buffer
            if char == b'\n':
                line = buffer.decode()[::-1]  # Reverse buffer to get line
                if "Score:" in line:
                    try:
                        score = float(line.split('Score: ')[1].split()[0])
                        scores.append(score)
                        lines_to_go -= 1
                    except ValueError:
                        # If there's an issue converting the score to a float, we ignore this line.
                        pass
                buffer = bytearray()
                if file.tell() == 1:
                    # If we are at the start of the file, read the last line.
                    file.seek(-1, os.SEEK_CUR)
            else:
                buffer.extend(char)
            file.seek(-1, os.SEEK_CUR)  # Move the pointer back to read the next char
        
        # Reverse the scores list so that it is in the correct order
        scores.reverse()
    return scores

# This function plots the scores in a bar chart.
def plot_scores_pf():
    # Define the number of groups and bar width
    n_groups = 20
    bar_width = 0.25

    # Create a position array for each set of bars
    index = np.arange(5)

    # Plot the bars
    fig, ax = plt.subplots(figsize=(15, 8))

    # Each set of bars is offset by the width of the bar
    bar1 = ax.bar(index, max_throughput_scores[5:10], bar_width, label='max_throughput_0')
    bar2 = ax.bar(index + bar_width, pf_scores[5:10], bar_width, label='pf_0')
    bar3 = ax.bar(index + bar_width * 2, mlwdf_scores[5:10], bar_width, label='mlwdf_0')

    # Add labels, title, and legend
    ax.set_xlabel('Slice Index')
    ax.set_ylabel('Scores')
    ax.set_title('Scores by slice - PF')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(range(1, n_groups+1))
    ax.legend()

    # Save the plot to a file
    plt.tight_layout()  # Adjust layout for better fit
    plt.savefig(f'./exp-customize-20slices/scores_comparison_pf_0.png')
    plt.close(fig)  # Close the plot to free up memory
    
    
def plot_scores_mt():
    # Define the number of groups and bar width
    n_groups = 20
    bar_width = 0.25

    # Create a position array for each set of bars
    index = np.arange(5)

    # Plot the bars
    fig, ax = plt.subplots(figsize=(15, 8))

    # Each set of bars is offset by the width of the bar
    bar1 = ax.bar(index, max_throughput_scores[:5], bar_width, label='max_throughput_0')
    bar2 = ax.bar(index + bar_width, pf_scores[:5], bar_width, label='pf_0')
    bar3 = ax.bar(index + bar_width * 2, mlwdf_scores[:5], bar_width, label='mlwdf_0')

    # Add labels, title, and legend
    ax.set_xlabel('Slice Index')
    ax.set_ylabel('Scores')
    ax.set_title('Scores by slice - MT')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(range(1, n_groups+1))
    ax.legend()

    # Save the plot to a file
    plt.tight_layout()  # Adjust layout for better fit
    plt.savefig(f'./exp-customize-20slices/scores_comparison_mt_0.png')
    plt.close(fig)  # Close the plot to free up memory
    
    
def plot_scores_mlwdf():
    # Define the number of groups and bar width
    n_groups = 20
    bar_width = 0.25

    # Create a position array for each set of bars
    index = np.arange(5)

    # Plot the bars
    fig, ax = plt.subplots(figsize=(15, 8))

    # Each set of bars is offset by the width of the bar
    bar1 = ax.bar(index, max_throughput_scores[10:15], bar_width, label='max_throughput_0')
    bar2 = ax.bar(index + bar_width, pf_scores[10:15], bar_width, label='pf_0')
    bar3 = ax.bar(index + bar_width * 2, mlwdf_scores[10:15], bar_width, label='mlwdf_0')

    # Add labels, title, and legend
    ax.set_xlabel('Slice Index')
    ax.set_ylabel('Scores')
    ax.set_title('Scores by slice - MLWD-F')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(range(1, n_groups+1))
    ax.legend()

    # Save the plot to a file
    plt.tight_layout()  # Adjust layout for better fit
    plt.savefig(f'./exp-customize-20slices/scores_comparison_mlwdf_0.png')
    plt.close(fig)  # Close the plot to free up memory


def plot_scores_video():
    # Define the number of groups and bar width
    n_groups = 20
    bar_width = 0.25

    # Create a position array for each set of bars
    index = np.arange(5)

    # Plot the bars
    fig, ax = plt.subplots(figsize=(15, 8))

    # Each set of bars is offset by the width of the bar
    bar1 = ax.bar(index, max_throughput_scores[10:15], bar_width, label='max_throughput_0')
    bar2 = ax.bar(index + bar_width, pf_scores[10:15], bar_width, label='pf_0')
    bar3 = ax.bar(index + bar_width * 2, mlwdf_scores[10:15], bar_width, label='mlwdf_0')

    # Add labels, title, and legend
    ax.set_xlabel('Slice Index')
    ax.set_ylabel('Scores')
    ax.set_title('Scores by slice - MLWD-F-VideoFlow')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(range(1, n_groups+1))
    ax.legend()

    # Save the plot to a file
    plt.tight_layout()  # Adjust layout for better fit
    plt.savefig(f'./exp-customize-20slices/scores_comparison_mlwdf_video_0.png')
    plt.close(fig)  # Close the plot to free up memory

# Process each log file and plot the scores
# for log_file in glob.glob(log_files_path):
#     file_name = os.path.basename(log_file)
#     scores = extract_last_20_scores(log_file)
#     plot_scores(scores, file_name)

# Plot the scores for the last 20 slices of each log file
max_throughput_scores = extract_last_20_scores('./exp-customize-20slices/max_throughput_0.log')
pf_scores = extract_last_20_scores('./exp-customize-20slices/pf_0.log')
mlwdf_scores = extract_last_20_scores('./exp-customize-20slices/mlwdf_0.log')
plot_scores_pf()
plot_scores_mt()
plot_scores_video()
print(max_throughput_scores)
print(pf_scores)
print(mlwdf_scores)


